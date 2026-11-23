"""LONAPEG-FAERS: per-quarter aggregation on retained (deduplicated) reports, window 2021Q3-2026Q3.
Reuses TESA-FAERS reduced quarterly parquet and retained_reports.parquet (same deduplication).
Writes partial files to data/processed/agg: report counts and PT counts by strata (all, US, HCP, pediatric <18, adult >=18),
prespecified-group and non-clinical-subcategory background counts, PS drug-PT pairs and drug marginals (MGPS),
and full extracts for reports naming any growth-hormone product (lonapegsomatropin + comparators)."""
import pandas as pd, numpy as np, sys, os
from lonaconfig import GROUPS, in_group, GH_ALL, nonclinical, nc_subcats
HERE=os.path.dirname(os.path.abspath(__file__)); TESA=f'{HERE}/../../tesamorelin-faers-2026-11/data/processed'
P=f'{TESA}/faers_reduced'; A=f'{HERE}/../data/processed/agg'; os.makedirs(A, exist_ok=True)
RET=f'{TESA}/retained_reports.parquet'
HCP={'MD','HP','PH','RN'}
gcache={}
def labels(pt):
    if pt not in gcache:
        g=[x for x in GROUPS if in_group(pt,x)]
        if nonclinical(pt): g += ['NC: any non-clinical'] + ['NC: '+s for s in nc_subcats(pt)]
        gcache[pt]=g
    return gcache[pt]
def num(df): df['primaryid']=pd.to_numeric(df.primaryid, errors='coerce'); return df
AGEF={'YR':1,'DEC':10,'MON':1/12,'WK':1/52,'DY':1/365,'HR':1/8760}
for q in sys.argv[1:]:
    if os.path.exists(f'{A}/{q}_talldrug.parquet'): print('skip',q); continue
    pids=set(pd.read_parquet(RET, columns=['primaryid'], filters=[('quarter','==',q)]).primaryid)
    demo=num(pd.read_parquet(f'{P}/{q}_demo.parquet')); demo=demo[demo.primaryid.isin(pids)].copy()
    reac=num(pd.read_parquet(f'{P}/{q}_reac.parquet')); reac['pt']=reac.pt.fillna('').str.strip().astype('category')
    reac=reac[reac.primaryid.isin(pids)&(reac.pt!='')]
    mixed=reac.loc[reac.pt!=reac.pt.str.upper(),'pt'].value_counts().rename('k').reset_index()
    reac['pt']=reac.pt.astype(str).str.upper().astype('category'); reac=reac.drop_duplicates(['primaryid','pt'])
    ctry=demo.occr_country.fillna(demo.reporter_country)
    age=pd.to_numeric(demo.age, errors='coerce')*demo.age_cod.map(AGEF)
    age=age.where((age>=0)&(age<=120))
    us=set(demo.loc[ctry.isin(['US','UNITED STATES']),'primaryid']); hcp=set(demo.loc[demo.occp_cod.isin(HCP),'primaryid'])
    ped=set(demo.loc[age<18,'primaryid']); adult=set(demo.loc[age>=18,'primaryid'])
    r_pids=set(reac.primaryid.unique())
    strata={'all':r_pids,'us':r_pids&us,'hcp':r_pids&hcp,'ped':r_pids&ped,'adult':r_pids&adult}
    stats=pd.DataFrame({'quarter':[q], **{f'N_{k}':[len(v)] for k,v in strata.items()}})
    pc=reac.groupby('pt',observed=True).primaryid.nunique().rename('n_all').to_frame()
    for k in ['us','hcp','ped','adult']: pc[f'n_{k}']=reac[reac.primaryid.isin(strata[k])].groupby('pt',observed=True).primaryid.nunique()
    pc=pc.fillna(0).astype('int64').reset_index()
    lab={p:labels(p) for p in reac.pt.cat.categories}
    hasg={p for p,v in lab.items() if v}
    gr=reac[reac.pt.isin(hasg)][['primaryid','pt']].copy(); gr['groups']=gr.pt.astype(str).map(lab)
    gr=gr[['primaryid','groups']].explode('groups').drop_duplicates()
    gbg=gr.groupby('groups').primaryid.apply(set)
    grows=[dict(quarter=q, group=g, **{f'n_{k}':len(s & v) for k,v in strata.items()}) for g,s in gbg.items()]
    gbg=pd.DataFrame(grows)
    drug=num(pd.read_parquet(f'{P}/{q}_drug.parquet', columns=['primaryid','role_cod','dkey'])); drug['dkey']=drug.dkey.astype('category'); drug['role_cod']=drug.role_cod.astype('category')
    dq=drug[drug.primaryid.isin(pids)]; del drug
    ps=dq[(dq.role_cod=='PS')&(dq.dkey!='')][['primaryid','dkey']].drop_duplicates(); ps['dkey']=ps.dkey.astype(str); ps=ps[ps.primaryid.isin(r_pids)]
    stats['N_ps_reports']=ps.primaryid.nunique()
    marg=ps.groupby('dkey').primaryid.nunique().rename('n_drug').reset_index()
    pairs=ps.merge(reac[['primaryid','pt']].assign(pt=reac.pt.astype(str)), on='primaryid').groupby(['dkey','pt']).size().rename('n').reset_index()
    del ps
    ghk=[k for k in dq.dkey.cat.categories if GH_ALL.search(k)]; tp=set(dq.loc[dq.dkey.isin(ghk), 'primaryid'])
    t=num(pd.read_parquet(f'{P}/{q}_tdrug.parquet')); t=t[t.primaryid.isin(tp)]
    outc=num(pd.read_parquet(f'{P}/{q}_outc.parquet')); ther=num(pd.read_parquet(f'{P}/{q}_tther.parquet')); indi=num(pd.read_parquet(f'{P}/{q}_tindi.parquet'))
    tdemo=demo[demo.primaryid.isin(tp)].copy()
    for name, df in [('ptdisplay',mixed),('stats',stats),('ptcounts',pc),('gbg',gbg),('marg',marg),('pairs',pairs),('tdrug',t),('tdemo',tdemo),
                     ('treac',reac[reac.primaryid.isin(tp)][['primaryid','pt']]),('toutc',outc[outc.primaryid.isin(tp)]),('tther',ther[ther.primaryid.isin(tp)]),
                     ('tindi',indi[indi.primaryid.isin(tp)]),('talldrug',dq[dq.primaryid.isin(tp)])]:
        df.astype({c:str for c in df.columns if df[c].dtype==object}).to_parquet(f'{A}/{q}_{name}.parquet', index=False)
    print(q, len(pids), len(r_pids), len(pairs), len(tp), flush=True)
