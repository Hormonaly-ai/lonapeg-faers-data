"""LONAPEG-FAERS (reused from TESA-FAERS): fit the DuMouchel (1999) two-component gamma-Poisson shrinkage prior to the PS drug x PT table.
Protocol section 8: zero-truncated likelihood first (non-zero cells only); if degenerate (alpha1 < 1e-3, or P < 0.001 or > 0.999)
or non-convergent, the non-truncated likelihood (zero cells included, as in DuMouchel's original MGPS) is used.
Data squashing: cells binned by rounded log10(n_drug) and log10(n_pt) (0.01); zero-cell counts obtained by
convolving the marginal histograms and subtracting observed cells."""
import pandas as pd, numpy as np, json, os, sys
from scipy.special import gammaln
from scipy.optimize import minimize
HERE=os.path.dirname(os.path.abspath(__file__)); D=f'{HERE}/../data/processed'
N=json.load(open(f'{D}/database_totals.json'))['N_ps_reports']
if not os.path.exists(f'{D}/mgps_pairs_with_E.parquet'):
    pairs=pd.read_parquet(f'{D}/mgps_pairs_ps.parquet'); marg=pd.read_parquet(f'{D}/mgps_drug_marginals.parquet')
    pc=pd.read_parquet(f'{D}/pt_counts_all_reports.parquet')[['pt','n_all']]
    pairs=pairs.merge(marg,on='dkey').merge(pc,on='pt'); pairs['E']=pairs.n_drug*pairs.n_all/N
    pairs.to_parquet(f'{D}/mgps_pairs_with_E.parquet',index=False)
pairs=pd.read_parquet(f'{D}/mgps_pairs_with_E.parquet',columns=['n','n_drug','n_all'])
marg=pd.read_parquet(f'{D}/mgps_drug_marginals.parquet'); pc=pd.read_parquet(f'{D}/pt_counts_all_reports.parquet')
R=100
bi=np.round(np.log10(marg.n_drug.values)*R).astype(int); bj=np.round(np.log10(pc.n_all.values[pc.n_all.values>0])*R).astype(int)
hi=np.bincount(bi); hj=np.bincount(bj); allc=np.convolve(hi,hj)  # index = bi+bj
pairs['b']=np.round(np.log10(pairs.n_drug)*R).astype(int)+np.round(np.log10(pairs.n_all)*R).astype(int)
obs_b=np.bincount(pairs.b.values,minlength=len(allc))
zero=allc-obs_b[:len(allc)]
assert (zero>=0).all()
lN=np.log10(N)
sq=pairs.groupby(['n','b']).size().rename('w').reset_index()
zb=np.nonzero(zero)[0]
n=np.r_[sq.n.values, np.zeros(len(zb))].astype(float)
E=np.r_[10**(sq.b.values/R-lN), 10**(zb/R-lN)]
w=np.r_[sq.w.values, zero[zb]].astype(float)
print('cells',w.sum(),'nonzero cells',sq.w.sum(),'squashed points',len(n),flush=True)
def lnb(n,a,b,E): return gammaln(a+n)-gammaln(a)-gammaln(n+1)+a*np.log(b/(b+E))+n*np.log(E/(b+E))
nz=n>0
def nll_trunc(t):
    a1,b1,a2,b2=np.exp(t[:4]); P=1/(1+np.exp(-t[4])); nn,EE,ww=n[nz],E[nz],w[nz]
    lf=np.logaddexp(np.log(P)+lnb(nn,a1,b1,EE),np.log1p(-P)+lnb(nn,a2,b2,EE))
    l0=np.logaddexp(np.log(P)+lnb(0,a1,b1,EE),np.log1p(-P)+lnb(0,a2,b2,EE))
    v=-np.sum(ww*(lf-np.log(-np.expm1(l0)))); return v if np.isfinite(v) else 1e300
trunc=[]
for s in [[0.2,0.1,2,4,1/3],[0.5,0.5,3,3,0.2],[1,1,5,5,0.5]]:
    t0=np.r_[np.log(s[:4]),np.log(s[4]/(1-s[4]))]
    r=minimize(nll_trunc,t0,method='Nelder-Mead',options=dict(maxiter=3000,xatol=1e-6,fatol=1e-3))
    trunc.append(dict(start=s,nll=float(r.fun),converged=bool(r.success),par=list(np.exp(r.x[:4]))+[float(1/(1+np.exp(-r.x[4])))])); print('trunc',trunc[-1],flush=True)
bt=min(trunc,key=lambda x:x['nll']); a1t,b1t,a2t,b2t,Pt=bt['par']
degenerate = (a1t<1e-3) or (Pt<1e-3) or (Pt>0.999) or (not bt['converged'])
def nll(t):
    a1,b1,a2,b2=np.exp(t[:4]); P=1/(1+np.exp(-t[4]))
    lf=np.logaddexp(np.log(P)+lnb(n,a1,b1,E),np.log1p(-P)+lnb(n,a2,b2,E))
    v=-np.sum(w*lf); return v if np.isfinite(v) else 1e300
starts=[[0.2,0.1,2,4,1/3],[0.5,0.5,3,3,0.2],[1,1,5,5,0.5],[0.05,0.05,1,1,0.1]]
best=None; fits=[]
for s in starts:
    t0=np.r_[np.log(s[:4]),np.log(s[4]/(1-s[4]))]
    r=minimize(nll,t0,method='Nelder-Mead',options=dict(maxiter=3000,xatol=1e-6,fatol=1e-3))
    r=minimize(nll,r.x,method='BFGS')
    fits.append(dict(start=s,nll=float(r.fun),par=list(np.exp(r.x[:4]))+[float(1/(1+np.exp(-r.x[4])))]))
    print(fits[-1],flush=True)
    if best is None or r.fun<best.fun: best=r
a1,b1,a2,b2=np.exp(best.x[:4]); P=1/(1+np.exp(-best.x[4]))
if not degenerate: a1,b1,a2,b2,P=a1t,b1t,a2t,b2t,Pt
prior=dict(truncated_fits=trunc, truncated_degenerate=bool(degenerate),alpha1=float(a1),beta1=float(b1),alpha2=float(a2),beta2=float(b2),P=float(P),negloglik=float(best.fun),
           n_cells_total=float(w.sum()),n_nonzero_cells=int(sq.w.sum()),n_squashed=int(len(n)),N=int(N),starts=fits,
           likelihood='non-truncated, zero cells included' if degenerate else 'zero-truncated')
json.dump(prior,open(f'{D}/mgps_prior.json','w'),indent=1); print({k:prior[k] for k in ['alpha1','beta1','alpha2','beta2','P']})
