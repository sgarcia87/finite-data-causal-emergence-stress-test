#!/usr/bin/env python3
"""
Final Emergence Stress Benchmark 1.0

Purpose
-------
Retain only findings that survived earlier falsification attempts:

1) Exact behavior-preserving state refinement can change rank-based clear-CE
   status while leaving lumped dynamics unchanged.
2) Exact empirical rank is unstable under finite transition sampling.
3) Forced cutoff/rank selectors can compress full-rank systems.
4) A one-sided "resolved modes" claim is safer than asserting exact rank or the
   absence of weaker modes.

This script does NOT claim to solve causal-emergence inference.
"""

import numpy as np
import pandas as pd
import math

SEED=2026091414
NCHAINS=150
NSPLITS=60
ROW_SAMPLES=[50,100,200,500,1000,5000]
N=16
rng=np.random.default_rng(SEED)

def EI(P):
    P=np.asarray(P,float)
    n=len(P)
    joint=P/n
    py=joint.sum(0)
    out=0.0
    for i in range(n):
        for j in range(n):
            q=joint[i,j]
            if q>0 and py[j]>0:
                out += q*np.log2(P[i,j]/py[j])
    return float(out)

def clone_equal(P,Nout=16):
    P=np.asarray(P,float)
    r=len(P)
    m=Nout//r
    assert r*m==Nout
    par=np.repeat(np.arange(r),m)
    Q=np.zeros((Nout,Nout),float)
    for a,i in enumerate(par):
        for b,j in enumerate(par):
            Q[a,b]=P[i,j]/m
    return Q,par

def lump(Q,par,r):
    L=np.zeros((r,r),float)
    for i in range(r):
        rows=np.where(par==i)[0]
        row=Q[rows].mean(0)
        for j in range(r):
            L[i,j]=row[par==j].sum()
    return L

def sample_counts(P,m):
    C=np.zeros_like(P,dtype=int)
    for i in range(len(P)):
        C[i]=rng.multinomial(m,P[i])
    return C

def empirical_tpm(C):
    s=C.sum(1,keepdims=True)
    X=np.zeros_like(C,dtype=float)
    good=s[:,0]>0
    X[good]=C[good]/s[good]
    return X

def linear_elbow_rank(s):
    y=np.asarray(s,float)
    x=np.linspace(0,1,len(y))
    yy=(y-y.min())/(y.max()-y.min()+1e-300)
    den=math.hypot(yy[-1]-yy[0],1.0)
    d=np.abs((yy[-1]-yy[0])*x-yy+yy[0])/(den+1e-300)
    return int(np.argmax(d[1:-1])+2)

def forced_rank(s,method):
    if method=="log_gap":
        g=np.log(np.maximum(s[:-1],1e-300))-np.log(np.maximum(s[1:],1e-300))
        return int(np.argmax(g)+1)
    if method=="linear_elbow":
        return linear_elbow_rank(s)
    if method=="svht":
        # Generic square-matrix hard-threshold baseline.
        # It is NOT presented as a method specifically validated for TPMs.
        eps=2.858*np.median(s)
        return int((s>eps).sum())
    raise ValueError(method)

def numerical_rank(s,reltol=1e-13):
    return int((s>reltol*s[0]).sum()) if s[0]>0 else 0

def resolved_lower_bound(C):
    """
    Cross-split resolved-mode lower bound.

    Per split:
      D=(P_A-P_B)/2
      G=(P_A^T P_B + P_B^T P_A)/2
      H=D^T D
      noise=lambda_max(H)
      margin_i=(lambda_i(G)-noise)/lambda_1(G)

    A leading mode is counted as resolved only when its 5th percentile margin
    across splits is >0.

    This returns ONLY a lower bound / resolved-mode count. It does not assert
    that unresolved modes are absent.
    """
    S,n=NSPLITS,C.shape[0]
    A=rng.binomial(C[None,:,:],0.5,size=(S,n,n))
    B=C[None,:,:]-A

    sa=A.sum(2,keepdims=True)
    sb=B.sum(2,keepdims=True)
    good=(sa[:,:,0]>0).all(1)&(sb[:,:,0]>0).all(1)
    if good.sum()<20:
        return np.nan

    A=A[good]/sa[good]
    B=B[good]/sb[good]

    D=(A-B)/2
    G=(np.matmul(np.transpose(A,(0,2,1)),B)+
       np.matmul(np.transpose(B,(0,2,1)),A))/2
    H=np.matmul(np.transpose(D,(0,2,1)),D)

    eigG=np.linalg.eigvalsh(G)[:,::-1]
    eigH=np.linalg.eigvalsh(H)[:,::-1]
    noise=np.maximum(eigH[:,0],0.0)
    lam1=np.maximum(eigG[:,0],1e-15)

    margins=(eigG-noise[:,None])/lam1[:,None]
    q05=np.quantile(margins,.05,axis=0)

    r=0
    for v in q05:
        if v>0:
            r+=1
        else:
            break
    return int(r)

def main():
    # Panel A: exact representation refinement
    rep_rows=[]
    for seed in range(NCHAINS):
        P=rng.dirichlet(np.ones(4)*1.5,size=4)
        Q,par=clone_equal(P,16)
        L=lump(Q,par,4)
        sP=np.linalg.svd(P,compute_uv=False)
        sQ=np.linalg.svd(Q,compute_uv=False)

        rep_rows.append({
            "seed":seed,
            "EI_base":EI(P),
            "EI_equal_clone":EI(Q),
            "EI_difference":EI(Q)-EI(P),
            "base_N":4,
            "clone_N":16,
            "base_rank":numerical_rank(sP),
            "clone_rank":numerical_rank(sQ),
            "base_clear_CE":int(numerical_rank(sP)<4),
            "clone_clear_CE":int(numerical_rank(sQ)<16),
            "lumping_error":float(np.abs(L-P).max())
        })
    rep=pd.DataFrame(rep_rows)

    # Panels B/C
    rows=[]

    def record_family(P,family,truth_rank,m):
        C=sample_counts(P,m)
        X=empirical_tpm(C)
        s=np.linalg.svd(X,compute_uv=False)
        rr=resolved_lower_bound(C)

        row={
            "family":family,
            "truth_rank":truth_rank,
            "samples_per_row":m,
            "total_transitions":N*m,
            "numerical_rank_empirical":numerical_rank(s),
            "clear_CE_empirical":int(numerical_rank(s)<N),
            "resolved_lower_bound":rr,
            "lower_bound_valid":int(rr<=truth_rank) if np.isfinite(rr) else 0,
            "resolved_fraction_of_truth":float(rr/truth_rank) if np.isfinite(rr) else np.nan,
        }
        for method in ["log_gap","linear_elbow","svht"]:
            k=forced_rank(s,method)
            row[f"{method}_rank"]=k
            row[f"{method}_exact"]=int(k==truth_rank)
            row[f"{method}_compresses"]=int(k<N)
        return row

    # Exact low-rank families
    for truth in [2,4,8]:
        for seed in range(NCHAINS):
            P0=rng.dirichlet(np.ones(truth)*1.5,size=truth)
            Q,_=clone_equal(P0,16)
            for m in ROW_SAMPLES:
                row=record_family(Q,f"EXACT_RANK_{truth}",truth,m)
                row["seed"]=seed
                rows.append(row)

    # Well-conditioned full-rank controls
    for alpha in [.2,.4]:
        P=(1-alpha)*np.eye(N)+alpha*np.ones((N,N))/N
        for seed in range(NCHAINS):
            for m in ROW_SAMPLES:
                row=record_family(P,f"FULL_RANK_WELL_a{alpha}",16,m)
                row["seed"]=seed
                rows.append(row)

    # Near-low-rank but mathematically full-rank
    for delta in [.001,.01,.05]:
        for seed in range(NCHAINS):
            P4=rng.dirichlet(np.ones(4)*1.5,size=4)
            Q,_=clone_equal(P4,16)
            R=rng.dirichlet(np.ones(N)*1.5,size=N)
            Xtrue=(1-delta)*Q+delta*R
            Xtrue/=Xtrue.sum(1,keepdims=True)
            truth=numerical_rank(np.linalg.svd(Xtrue,compute_uv=False))

            for m in ROW_SAMPLES:
                row=record_family(Xtrue,f"NEAR_LOW_FULL_d{delta}",truth,m)
                row["seed"]=seed
                rows.append(row)

    bench=pd.DataFrame(rows)

    low=bench[bench.family.str.startswith("EXACT_RANK_")].copy()
    low_summary=(low.groupby(["family","samples_per_row"],as_index=False)
                 .agg(
                     clear_CE_detection_rate=("clear_CE_empirical","mean"),
                     log_gap_exact_rate=("log_gap_exact","mean"),
                     linear_elbow_exact_rate=("linear_elbow_exact","mean"),
                     svht_exact_rate=("svht_exact","mean"),
                     lower_bound_valid_rate=("lower_bound_valid","mean"),
                     median_resolved_fraction=("resolved_fraction_of_truth","median"),
                     p10_resolved_fraction=("resolved_fraction_of_truth",lambda x:x.quantile(.1)),
                     p90_resolved_fraction=("resolved_fraction_of_truth",lambda x:x.quantile(.9)),
                 ))

    full=bench[~bench.family.str.startswith("EXACT_RANK_")].copy()
    full_summary=(full.groupby(["family","samples_per_row"],as_index=False)
                  .agg(
                      clear_CE_false_positive_rate=("clear_CE_empirical","mean"),
                      log_gap_false_compression_rate=("log_gap_compresses","mean"),
                      linear_elbow_false_compression_rate=("linear_elbow_compresses","mean"),
                      svht_false_compression_rate=("svht_compresses","mean"),
                      lower_bound_valid_rate=("lower_bound_valid","mean"),
                      median_resolved_lower_bound=("resolved_lower_bound","median"),
                      p90_resolved_lower_bound=("resolved_lower_bound",lambda x:x.quantile(.9)),
                  ))

    pub=[]
    for m in [100,500,5000]:
        ql=low[low.samples_per_row==m]
        qf=full[full.samples_per_row==m]
        pub.append({
            "samples_per_row":m,
            "lowrank_clear_CE_detection":ql.clear_CE_empirical.mean(),
            "lowrank_loggap_exact":ql.log_gap_exact.mean(),
            "lowrank_elbow_exact":ql.linear_elbow_exact.mean(),
            "lowrank_svht_exact":ql.svht_exact.mean(),
            "lower_bound_valid_all_families":bench[bench.samples_per_row==m].lower_bound_valid.mean(),
            "lowrank_median_resolved_fraction":ql.resolved_fraction_of_truth.median(),
            "fullrank_loggap_false_compression":qf.log_gap_compresses.mean(),
            "fullrank_elbow_false_compression":qf.linear_elbow_compresses.mean(),
            "fullrank_svht_false_compression":qf.svht_compresses.mean(),
            "fullrank_median_resolved_lower_bound":qf.resolved_lower_bound.median(),
        })
    pub=pd.DataFrame(pub)

    rep.to_csv("final_benchmark_representation_invariance.csv",index=False)
    bench.to_csv("final_benchmark_all_results.csv",index=False)
    low_summary.to_csv("final_benchmark_lowrank_summary.csv",index=False)
    full_summary.to_csv("final_benchmark_fullrank_summary.csv",index=False)
    pub.to_csv("final_benchmark_publication_endpoints.csv",index=False)

    print("=== PANEL A ===")
    print(pd.DataFrame([{
        "n_chains":len(rep),
        "max_lumping_error":rep.lumping_error.max(),
        "median_abs_EI_change":rep.EI_difference.abs().median(),
        "fraction_base_clear_CE":rep.base_clear_CE.mean(),
        "fraction_clone_clear_CE":rep.clone_clear_CE.mean(),
        "base_rank_median":rep.base_rank.median(),
        "clone_rank_median":rep.clone_rank.median(),
    }]).to_string(index=False))

    print("\n=== PUBLICATION ENDPOINTS ===")
    print(pub.to_string(index=False))

if __name__=="__main__":
    main()
