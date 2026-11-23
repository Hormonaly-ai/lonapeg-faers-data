"""LONAPEG-FAERS (reused from TESA-FAERS tesastats.py): disproportionality statistics (ROR, PRR, BCPNN IC, MGPS EBGM/EB05)."""
import numpy as np, json, os
from scipy.special import gammaln, digamma
from scipy.stats import gamma, chi2_contingency
from scipy.optimize import brentq
HERE=os.path.dirname(os.path.abspath(__file__))
PRIOR=json.load(open(f'{HERE}/../data/processed/mgps_prior.json'))
def _lnb(n,a,b,E): return gammaln(a+n)-gammaln(a)-gammaln(n+1)+a*np.log(b/(b+E))+n*np.log(E/(b+E))
def ebgm(n,E,p=PRIOR):
    a1,b1,a2,b2,P=p['alpha1'],p['beta1'],p['alpha2'],p['beta2'],p['P']
    l1=np.log(P)+_lnb(n,a1,b1,E); l2=np.log1p(-P)+_lnb(n,a2,b2,E); Q=1/(1+np.exp(l2-l1))
    elog=Q*(digamma(a1+n)-np.log(b1+E))+(1-Q)*(digamma(a2+n)-np.log(b2+E))
    cdf=lambda x: Q*gamma.cdf(x,a1+n,scale=1/(b1+E))+(1-Q)*gamma.cdf(x,a2+n,scale=1/(b2+E))-0.05
    eb05=brentq(cdf,1e-8,1e5)
    return float(np.exp(elog)), float(eb05)
def table2x2(a,b,c,d):
    out={}
    aa,bb,cc,dd=[x+0.5 for x in (a,b,c,d)] if min(a,b,c,d)==0 else (a,b,c,d)
    ror=(aa/bb)/(cc/dd); se=np.sqrt(1/aa+1/bb+1/cc+1/dd)
    out.update(ROR=ror,ROR_lo=np.exp(np.log(ror)-1.96*se),ROR_hi=np.exp(np.log(ror)+1.96*se))
    prr=(aa/(aa+bb))/(cc/(cc+dd)); sep=np.sqrt(1/aa-1/(aa+bb)+1/cc-1/(cc+dd))
    out.update(PRR=prr,PRR_lo=np.exp(np.log(prr)-1.96*sep),PRR_hi=np.exp(np.log(prr)+1.96*sep))
    try: out['chi2_yates']=float(chi2_contingency([[a,b],[c,d]],correction=True)[0])
    except ValueError: out['chi2_yates']=np.nan
    N=a+b+c+d; E=(a+b)*(a+c)/N
    out['E']=E; out['IC']=np.log2((a+0.5)/(E+0.5)); out['IC025']=out['IC']-3.3*(a+0.5)**-0.5-2*(a+0.5)**-1.5
    return out
def vs_all(a,nT,nPT,N):
    """a: target reports with PT; nT: target reports; nPT: all reports with PT; N: all reports."""
    b=nT-a; c=nPT-a; d=N-nT-c
    o=table2x2(a,b,c,d); o['EBGM'],o['EB05']=ebgm(a,nT*nPT/N); o['E_mgps']=nT*nPT/N
    o['sig_ROR']=bool(a>=3 and o['ROR_lo']>1); o['sig_PRR']=bool(a>=3 and o['PRR']>=2 and o['chi2_yates']>=4)
    o['sig_IC']=bool(o['IC025']>0); o['sig_EB05']=bool(o['EB05']>=2)
    o['signal_primary']=o['sig_ROR'] and o['sig_IC']; o['signal_highspec']=o['signal_primary'] and o['sig_EB05']
    return o
