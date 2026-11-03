#!/usr/bin/env python3
"""Bounded numerical checks for one Gaussian source and a five-point target.

Gauss--Hermite quadrature approximates the continuous Gaussian integrals.
Two quadrature resolutions are compared; this is not a proof or a claim that
the finite quadrature source itself satisfies a Euclidean Poincare inequality.
"""
import argparse
import json
import math
import platform
import time
from pathlib import Path

import numpy as np
import scipy
from scipy.linalg import eigh, null_space
from scipy.optimize import root
from scipy.special import roots_hermitenorm, logsumexp

Y=np.array([-1.0,-0.2,0.3,0.8,1.0])
NU=np.array([0.1,0.15,0.2,0.25,0.3])
LOGNU=np.log(NU)
B=null_space(np.ones((1,len(Y))))
POTENTIALS=[
    np.zeros(5),
    np.array([-0.3,0.2,0.1,-0.2,0.4]),
    np.array([0.7,-0.4,0.3,-0.6,0.2]),
    np.array([-1.0,-0.5,0.0,0.5,1.0]),
    np.array([0.2,0.2,-0.4,0.2,0.2]),
]


class GaussianModel:
    def __init__(self,n,lam):
        x,w=roots_hermitenorm(n)
        self.x=x
        self.w=w/math.sqrt(2*math.pi)
        self.cost=np.sin(x[:,None])*Y[None,:]
        self.lam=lam

    def evaluate(self,v):
        z=LOGNU[None,:]+(v[None,:]-self.cost)/self.lam
        lz=logsumexp(z,axis=1)
        P=np.exp(z-lz[:,None])
        q=self.w@P
        J=self.lam*(self.w@lz)-NU@v
        cross=P.T@(self.w[:,None]*P)
        D=np.diag(q)-cross
        C=cross-np.outer(q,q)
        return float(J),q,D,C

    def update(self,v):
        _,q,_,_=self.evaluate(v)
        vn=v+self.lam*np.log(NU/q)
        return vn-NU@vn

    def optimum(self):
        # Solve four independent marginal equations with the fifth potential fixed.
        def fun(z):
            return self.evaluate(np.r_[z,0.0])[1][:-1]-NU[:-1]
        def jac(z):
            return self.evaluate(np.r_[z,0.0])[2][:-1,:-1]/self.lam
        sol=root(fun,np.zeros(4),jac=jac,tol=1e-11)
        v=np.r_[sol.x,0.0]
        v-=NU@v
        J,q,_,_=self.evaluate(v)
        residual=float(np.max(np.abs(q-NU)))
        if residual>1e-11:
            raise RuntimeError(f"Optimum residual {residual} at lambda={self.lam}")
        # Recheck the exact discrete Sinkhorn update at this same quadrature optimum.
        stationary_error=float(np.max(np.abs(self.update(v)-v)))
        return v,J,residual,stationary_error,bool(sol.success)


def main(outdir,resolutions):
    started=time.monotonic()
    outdir.mkdir(parents=True,exist_ok=True)
    rep={"scope":"Floating-point Gauss--Hermite approximation, not a proof.",
         "source":"N(0,1)","cost":"sin(x)*y", "C_P":1.0,"L":1.0,
         "target_points":Y.tolist(),"target_weights":NU.tolist(),
         "resolutions":resolutions,"lambda_values":[0.5,0.2,0.1],
         "versions":{"python":platform.python_version(),"numpy":np.__version__,
                     "scipy":scipy.__version__},"runs":[],"comparison":[],"failures":[]}
    for n in resolutions:
        for lam in rep["lambda_values"]:
            model=GaussianModel(n,lam)
            vstar,Jstar,res,stationary,solver_success=model.optimum()
            gamma=lam*lam/(1+lam*lam)
            entries=[]
            vectors=POTENTIALS+[vstar+0.05*np.array([-0.4,0.7,-0.8,0.1,0.5])]
            for k,v0 in enumerate(vectors):
                v=v0-NU@v0
                J,q,D,C=model.evaluate(v)
                delta=J-Jstar
                kl=float(NU@np.log(NU/q))
                vn=model.update(v)
                nextdelta=model.evaluate(vn)[0]-Jstar
                error_ratio=delta/((lam+1/lam)*kl)
                contraction_ratio=nextdelta/delta
                # Poincare: C <= lambda^{-2} D. Check all target directions,
                # excluding the constant direction, by a generalised eigensolve.
                Db=B.T@D@B
                Cb=B.T@C@B
                variance_ratio=float(eigh(Cb,Db,eigvals_only=True)[-1])*lam*lam
                M=(D-gamma*(np.diag(q)-np.outer(q,q)))/lam
                min_hessian_eigenvalue=float(np.linalg.eigvalsh(B.T@M@B)[0])
                descent_ratio=lam*kl/(delta-nextdelta)
                entry={"potential_index":k,"potential":v.tolist(),"dual_gap":delta,
                       "marginal_kl":kl,"global_error_bound_ratio":error_ratio,
                       "next_dual_gap":nextdelta,"gap_contraction_ratio":contraction_ratio,
                       "theoretical_contraction_factor":1-gamma,
                       "conditional_variance_bound_ratio":variance_ratio,
                       "relative_hessian_min_eigenvalue":min_hessian_eigenvalue,
                       "descent_bound_ratio":descent_ratio}
                entries.append(entry)
                if delta<=0 or error_ratio>1+1e-7 or variance_ratio>1+1e-7 or min_hessian_eigenvalue < -1e-10:
                    rep["failures"].append([n,lam,k,"error/hessian inequality",entry])
                if contraction_ratio>1-gamma+1e-7 or descent_ratio>1+1e-7:
                    rep["failures"].append([n,lam,k,"contraction/descent",entry])
            # One actual nonzero-init trajectory, until a marginal residual target.
            v=POTENTIALS[1].copy();v-=NU@v
            trajectory=[]
            max_trajectory_ratio=0.0
            for step in range(20001):
                J,q,_,_=model.evaluate(v)
                delta=J-Jstar
                residual=float(np.max(np.abs(q-NU)))
                if step in [0,1,2,5,10,20,50,100,200,500,1000,2000,5000,10000,20000]:
                    trajectory.append({"step":step,"gap":delta,"marginal_residual":residual})
                if residual<1e-12:
                    break
                vn=model.update(v)
                nextdelta=model.evaluate(vn)[0]-Jstar
                # Restrict ratios to gaps above floating-point cancellation.
                if delta>1e-9:
                    ratio=nextdelta/delta
                    max_trajectory_ratio=max(max_trajectory_ratio,ratio)
                    if ratio>1-gamma+1e-7:
                        rep["failures"].append([n,lam,step,"trajectory contraction",ratio])
                v=vn
            rep["runs"].append({"nodes":n,"lambda":lam,
                "quadrature_normalisation":float(model.w.sum()),
                "optimum_potential":vstar.tolist(),"optimum_objective":Jstar,
                "optimum_marginal_residual":res,"optimum_update_error":stationary,
                "root_solver_success_flag":solver_success,"cases":entries,
                "trajectory_steps":step,"final_marginal_residual":residual,
                "max_trajectory_gap_ratio_above_1e-9":max_trajectory_ratio,
                "trajectory":trajectory})
    lo,hi=resolutions
    for lam in rep["lambda_values"]:
        r1=next(z for z in rep["runs"] if z["nodes"]==lo and z["lambda"]==lam)
        r2=next(z for z in rep["runs"] if z["nodes"]==hi and z["lambda"]==lam)
        fields=["dual_gap","marginal_kl","next_dual_gap",
                "conditional_variance_bound_ratio","relative_hessian_min_eigenvalue"]
        differences={f:max(abs(a[f]-b[f]) for a,b in zip(r1["cases"],r2["cases"])) for f in fields}
        differences["optimum_potential"]=float(np.max(np.abs(np.array(r1["optimum_potential"])-np.array(r2["optimum_potential"]))))
        rep["comparison"].append({"lambda":lam,"maximum_absolute_differences":differences})
        if max(differences.values())>2e-7:
            rep["failures"].append([lam,"quadrature resolution agreement",differences])
    rep["elapsed_seconds"]=time.monotonic()-started
    rep["passed"]=not rep["failures"]
    (outdir/"multitarget_gaussian_checks.json").write_text(json.dumps(rep,indent=2)+"\n")
    lines=[rep["scope"],f"Source N(0,1), c(x,y)=sin(x)y, five target atoms; C_P=L=1.",
           f"Gauss--Hermite resolutions: {resolutions}; lambda values: {rep['lambda_values']}."]
    for run in rep["runs"]:
        lines.append(f"n={run['nodes']}, lambda={run['lambda']}: optimum marginal residual "
            f"{run['optimum_marginal_residual']:.3g}; {len(run['cases'])} potential cases; "
            f"max error-bound ratio {max(z['global_error_bound_ratio'] for z in run['cases']):.9g}; "
            f"max variance-bound ratio {max(z['conditional_variance_bound_ratio'] for z in run['cases']):.9g}; "
            f"trajectory {run['trajectory_steps']} steps to residual {run['final_marginal_residual']:.3g}.")
    for comp in rep["comparison"]:
        lines.append(f"lambda={comp['lambda']}: resolution differences {comp['maximum_absolute_differences']}")
    lines.append(f"Failures: {len(rep['failures'])}; elapsed {rep['elapsed_seconds']:.3f}s.")
    for failure in rep["failures"]:
        lines.append("FAIL "+repr(failure))
    (outdir/"multitarget_gaussian_checks.log").write_text("\n".join(lines)+"\n")
    print("\n".join(lines))
    return 0 if rep["passed"] else 1


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--output-dir",type=Path,default=Path(__file__).resolve().parents[1]/"checks")
    parser.add_argument("--resolutions",type=int,nargs=2,default=[4096,16384])
    args=parser.parse_args()
    raise SystemExit(main(args.output_dir,args.resolutions))
