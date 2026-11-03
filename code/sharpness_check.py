#!/usr/bin/env python3
"""Independent floating-point checks for the two-point Sinkhorn examples.

These computations are numerical sanity checks, not certificates of the
unrestricted theorem. The analytical derivations are in
checks/independent_sharpness_audit.md in the bundle. Run with Python 3,
NumPy and SciPy; Matplotlib is optional for regenerating the plot.
"""
import argparse
import json
import math
import platform
import time
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import quad
from scipy.special import expit

LOG2 = math.log(2.0)
LAMBDAS = [0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001]
AMPLITUDES = [1e-5, 0.001, 0.01, 0.05, 0.2, 0.5, 0.8, 1.0, 2.0, 5.0]


def logcosh(z):
    a = abs(z)
    return a + math.log1p(math.exp(-2.0 * a)) - LOG2


def logsinh(z):
    a = abs(z)
    if a == 0:
        return -math.inf
    if a < 1.0:
        return math.log(math.sinh(a))
    return a + math.log1p(-math.exp(-2.0 * a)) - LOG2


def logsech2(z):
    return -2.0 * logcosh(z)


def logexpit(z):
    return -float(np.logaddexp(0.0, -z))


def half_density(source, x):
    """Density of |X| on (0,1), normalised to integrate to one."""
    if x <= 0 or x > 1.0 + 1e-12:
        return 0.0 if source == "log_cusp" else (1.0 if x == 0 else 0.0)
    x = min(x, 1.0)
    if source == "uniform":
        return 1.0
    return 0.8 * x * math.log(math.e / x) ** 2


def cp_bound(source):
    return 4.0 / math.pi**2 if source == "uniform" else 0.75


def segmented_quad(fun, lo, hi, points=()):
    if hi <= lo:
        return 0.0
    cuts = sorted(set([lo, hi] + [float(z) for z in points if lo < z < hi]))
    return sum(quad(fun, a, b, epsabs=2e-12, epsrel=2e-11, limit=150)[0]
               for a, b in zip(cuts, cuts[1:]))


def marginal_mean(source, lam, t):
    if t == 0:
        return 0.0
    sign = 1.0 if t > 0 else -1.0
    s = abs(t) / lam
    ln = logsinh(2.0 * s)
    lc = logcosh(2.0 * s)
    upper = min(1.0 / lam, max(40.0, s + 40.0))
    scale = s if s < 1.0 else 1.0

    def fun(u):
        pair = math.exp(ln - float(np.logaddexp(lc, logcosh(2.0 * u))))
        return lam * half_density(source, lam * u) * pair / scale

    return sign * scale * segmented_quad(fun, 0.0, upper, [1.0, s - 5.0, s, s + 5.0])


def marginal_logs(source, lam, t):
    """log(q(+1)), log(q(-1)), with endpoint scaling for rare targets."""
    if abs(t) <= 0.9:
        m = marginal_mean(source, lam, t)
        return math.log1p(m) - LOG2, math.log1p(-m) - LOG2
    at = abs(t)
    peak = logexpit(2.0 * (1.0 - at) / lam)
    upper = min(1.0 / lam, max(0.0, (1.0 - at) / lam) + 40.0)

    def fun(v):
        x = max(0.0, 1.0 - lam * v)
        left = logexpit(-2.0 * ((1.0 + at) / lam - v))
        right = logexpit(2.0 * ((1.0 - at) / lam - v))
        return 0.5 * lam * half_density(source, x) * (
            math.exp(left - peak) + math.exp(right - peak))

    scaled = segmented_quad(fun, 0.0, upper, [1.0, (1.0-at)/lam])
    logminor = peak + math.log(scaled)
    logmajor = math.log1p(-math.exp(logminor))
    return (logmajor, logminor) if t > 0 else (logminor, logmajor)


def marginal_kl(source, lam, t):
    if abs(t) <= 0.9:
        m = marginal_mean(source, lam, t)
        return -0.5 * math.log1p(-m * m)
    lp, lm = marginal_logs(source, lam, t)
    return -LOG2 - 0.5 * (lp + lm)


def gap(source, lam, t):
    if t == 0:
        return 0.0
    s = abs(t) / lam
    ls = 2.0 * logsinh(s)
    upper = min(1.0 / lam, max(40.0, s + 40.0))
    scale = s*s if s < 1.0 else 1.0

    def fun(u):
        paired = float(np.logaddexp(0.0, ls - 2.0 * logcosh(u)))
        return 0.5 * lam * lam * half_density(source, lam * u) * paired / scale

    return scale * segmented_quad(fun, 0.0, upper, [1.0, s - 5.0, s, s + 5.0])


def log_conditional_variance(source, lam, t):
    """Log E sech²((X+t)/lambda), resolving its narrow moving peak."""
    t = abs(t)
    lo = (t - 1.0) / lam
    hi = (t + 1.0) / lam
    peak_u = max(0.0, lo)
    peak = logsech2(peak_u)
    if lo > 0:
        hi = min(hi, lo + 40.0)
    else:
        lo, hi = max(lo, -40.0), min(hi, 40.0)

    def fun(u):
        x = min(1.0, abs(lam * u - t))
        return 0.5 * lam * half_density(source, x) * math.exp(logsech2(u)-peak)

    scaled = segmented_quad(fun, lo, hi, [-1.0, 0.0, 1.0, t/lam])
    return peak + math.log(scaled)


def update(source, lam, t):
    if abs(t) <= 0.9:
        return t - lam * math.atanh(marginal_mean(source, lam, t))
    lp, lm = marginal_logs(source, lam, t)
    return t - 0.5 * lam * (lp-lm)


def run(outdir):
    started = time.monotonic()
    outdir.mkdir(parents=True, exist_ok=True)
    report = {
        "scope": "Floating-point sanity checks only; not a proof certificate.",
        "versions": {"python": platform.python_version(), "numpy": np.__version__,
                     "scipy": scipy.__version__},
        "lambda_values": LAMBDAS,
        "potential_amplitudes": AMPLITUDES,
        "sources": {},
        "failures": [],
    }
    lines = [report["scope"]]
    for source in ["uniform", "log_cusp"]:
        norm = quad(lambda x: half_density(source, x), 0.0, 1.0,
                    epsabs=1e-12, epsrel=1e-12)[0]
        cp = cp_bound(source)
        records, detailed = [], []
        max_error_ratio = 0.0
        max_variance_ratio = 0.0
        max_descent_ratio = 0.0
        for lam in LAMBDAS:
            qgap = math.exp(log_conditional_variance(source, lam, 0.0))
            r = 1.0 - qgap
            rec = {"lambda": lam, "q_gap": qgap, "local_potential_factor": r,
                   "local_dual_gap_factor": r*r,
                   "poincare_gamma": lam*lam/(lam*lam+cp)}
            if source == "uniform":
                exact = lam * math.tanh(1.0/lam)
                rec["uniform_exact_q_gap"] = exact
                rec["uniform_quadrature_absolute_error"] = abs(qgap-exact)
                if abs(qgap-exact) > 2e-10:
                    report["failures"].append([source,lam,"uniform qgap"])
            else:
                lead = 0.8*LOG2*lam*lam*math.log(1.0/lam)**2
                rec["q_gap_over_leading_asymptotic"] = qgap/lead
            records.append(rec)
            A = cp/(lam*lam)
            C = lam+cp/lam
            for t in AMPLITUDES:
                delta = gap(source,lam,t)
                kl = marginal_kl(source,lam,t)
                ratio = delta/(C*kl)
                lp,lm = marginal_logs(source,lam,t)
                logtotal = math.log(4.0)+lp+lm
                logD = log_conditional_variance(source,lam,t)
                var_ratio = math.expm1(logtotal-logD)/A
                nxt=update(source,lam,t)
                nextgap=gap(source,lam,nxt)
                # Sinkhorn improvement >= lambda KL(nu || q).
                improvement=delta-nextgap
                descent_ratio=lam*kl/improvement if improvement > 0 else math.inf
                detailed.append({"lambda":lam,"t":t,"dual_gap":delta,
                                 "marginal_kl":kl,"error_bound_ratio":ratio,
                                 "conditional_variance_bound_ratio":var_ratio,
                                 "update":nxt,"next_gap":nextgap,
                                 "descent_bound_ratio":descent_ratio})
                max_error_ratio=max(max_error_ratio,ratio)
                max_variance_ratio=max(max_variance_ratio,var_ratio)
                max_descent_ratio=max(max_descent_ratio,descent_ratio)
                if ratio>1.0+1e-7:
                    report["failures"].append([source,lam,t,"global error bound",ratio])
                if var_ratio>1.0+1e-7 or var_ratio < -1e-6:
                    report["failures"].append([source,lam,t,"conditional variance",var_ratio])
                if descent_ratio>1.0+2e-6:
                    report["failures"].append([source,lam,t,"descent",descent_ratio])
                if not 0 < nxt < t:
                    report["failures"].append([source,lam,t,"monotone scalar update",nxt])

        trajectories=[]
        for lam in [0.5,0.2,0.1,0.05,0.02]:
            t=0.8
            samples=[]
            for n in range(15001):
                if n in [0,1,2,5,10,20,50,100,200,500,1000,2000,5000,10000,15000]:
                    samples.append({"iteration":n,"t":t,"gap":gap(source,lam,t)})
                if t <= lam*1e-4:
                    break
                t=update(source,lam,t)
            nxt=update(source,lam,t)
            observed=gap(source,lam,nxt)/gap(source,lam,t)
            r=1.0-math.exp(log_conditional_variance(source,lam,0.0))
            trajectories.append({"lambda":lam,"initial_t":0.8,"steps":n,
                                 "final_t":t,"observed_gap_ratio":observed,
                                 "predicted_gap_ratio":r*r,
                                 "ratio_absolute_error":abs(observed-r*r),
                                 "samples":samples})
            if abs(observed-r*r)>2e-7:
                report["failures"].append([source,lam,"trajectory ratio",observed,r*r])
        report["sources"][source]={"normalisation":norm,"poincare_constant_used":cp,
            "local_rates":records,"grid_checks":detailed,"trajectories":trajectories,
            "max_global_error_bound_ratio":max_error_ratio,
            "max_conditional_variance_bound_ratio":max_variance_ratio,
            "max_descent_bound_ratio":max_descent_ratio}
        lines.append(f"{source}: normalisation={norm:.16g}; C_P={cp:.16g}; "
                     f"max error-bound ratio={max_error_ratio:.9g}; "
                     f"max conditional-variance ratio={max_variance_ratio:.9g}; "
                     f"max descent-bound ratio={max_descent_ratio:.9g}.")
        lines.append(f"{source}: {len(detailed)} parameter-grid checks and "
                     f"{len(trajectories)} actual nonzero-initialisation trajectories.")
        for tr in trajectories:
            lines.append(f"  lambda={tr['lambda']}: {tr['steps']} steps; "
                         f"final ratio={tr['observed_gap_ratio']:.12g}; "
                         f"predicted={tr['predicted_gap_ratio']:.12g}; "
                         f"absolute difference={tr['ratio_absolute_error']:.3g}.")
    report["elapsed_seconds"]=time.monotonic()-started
    report["passed"] = not report["failures"]
    lines.append(f"Failures: {len(report['failures'])}. Elapsed: {report['elapsed_seconds']:.3f}s.")
    for failure in report["failures"]:
        lines.append("FAIL "+repr(failure))
    (outdir/"sharpness_checks.json").write_text(json.dumps(report,indent=2)+"\n")
    (outdir/"sharpness_checks.log").write_text("\n".join(lines)+"\n")
    print("\n".join(lines))

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig,ax=plt.subplots(figsize=(7.4,4.6),layout="constrained")
        fig.patch.set_facecolor("#faf7f2")
        ax.set_facecolor("#faf7f2")
        for name,label,colour in [("uniform","Uniform source (exact)","#b06535"),
                                  ("log_cusp","Fixed logarithmic cusp","#6d3d73")]:
            rows=report["sources"][name]["local_rates"]
            xx=np.array([r["lambda"] for r in rows])
            yy=np.array([1-r["local_dual_gap_factor"] for r in rows])
            ax.loglog(xx,yy,"o-",color=colour,label=label,lw=2,ms=5)
        xx=np.geomspace(min(LAMBDAS),max(LAMBDAS),160)
        ax.loglog(xx,xx**2/(xx**2+0.75),"--",color="#386f68",lw=1.8,
                  label="Poincaré guarantee, $C_P=3/4$")
        ax.set_xlabel("Entropy regularisation $\\lambda$")
        ax.set_ylabel("Asymptotic dual-gap improvement $1-r_\\lambda^2$")
        ax.set_title("A fixed source can approach the quadratic rate exponent")
        ax.grid(True,which="major",alpha=.16)
        ax.legend(frameon=False,loc="upper left")
        for spine in ax.spines.values():
            spine.set_color("#b8ada1")
        fig.savefig(outdir/"sharpness_plot.png",dpi=190)
        fig.savefig(outdir/"sharpness_plot.svg")
        plt.close(fig)
    except ImportError:
        print("Plot skipped: matplotlib unavailable.")
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--output-dir",type=Path,default=Path(__file__).resolve().parents[1]/"checks")
    args=parser.parse_args()
    raise SystemExit(run(args.output_dir))
