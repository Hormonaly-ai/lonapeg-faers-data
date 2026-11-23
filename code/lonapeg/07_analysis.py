"""LONAPEG-FAERS main analysis (protocol v1.0, OSF 10.17605/OSF.IO/UB6W2). Adapted from TESA-FAERS 07_analysis.py.
All manuscript numbers come from outputs/ written here."""
import pandas as pd, numpy as np, json, os, re
from lonaconfig import LONA, GROUPS, in_group, nonclinical, nc_subcats, NC_SUB, labeled, indication_related, drug_class
from lonastats import vs_all, table2x2
from scipy.stats import weibull_min, norm
HERE = os.path.dirname(os.path.abspath(__file__)); D = f'{HERE}/../data/processed'; T = f'{HERE}/../outputs/tables'
os.makedirs(T, exist_ok=True)
R = {}
tdrug = pd.read_parquet(f'{D}/target_tdrug.parquet'); tdemo = pd.read_parquet(f'{D}/target_tdemo.parquet')
treac = pd.read_parquet(f'{D}/target_treac.parquet'); toutc = pd.read_parquet(f'{D}/target_toutc.parquet')
tther = pd.read_parquet(f'{D}/target_tther.parquet'); tindi = pd.read_parquet(f'{D}/target_tindi.parquet')
talld = pd.read_parquet(f'{D}/target_talldrug.parquet')
for df in (tdrug, tdemo, treac, toutc, tther, tindi, talld): df['primaryid'] = pd.to_numeric(df.primaryid).astype('int64')
pc = pd.read_parquet(f'{D}/pt_counts_all_reports.parquet').set_index('pt')
gb = pd.read_csv(f'{D}/group_background_counts.csv').set_index('group')
DISP = pd.read_csv(f'{D}/pt_display_names.csv').set_index('PT').display.to_dict()
def disp(pt): return DISP.get(pt, pt.capitalize())
tot = json.load(open(f'{D}/database_totals.json')); N = tot['N_reports_with_pt']
R['database'] = dict(tot, dedup_source='TESA-FAERS retained_reports.parquet (latest CASEVERSION per CASEID across 2010Q4-2026Q3, deleted cases removed)')
R['database']['tesa_dedup_summary'] = json.load(open(f'{HERE}/../../tesamorelin-faers-2026-11/data/processed/dedup_summary.json'))
R['mgps_prior'] = {k: v for k, v in json.load(open(f'{D}/mgps_prior.json')).items() if k in ['alpha1','beta1','alpha2','beta2','P','likelihood','truncated_degenerate']}
withpt = set(treac.primaryid)
# ---------- case definition ----------
talld['cls'] = talld.dkey.map(drug_class)
tdrug['is_L'] = tdrug.drugname.str.contains(LONA) | tdrug.prod_ai.str.contains(LONA)
lname = tdrug[tdrug.is_L][['primaryid','role_cod']]
# Daily somatropin vs Serostim/Zorbtive must use DRUGNAME (PROD_AI is 'SOMATROPIN' for all); tdrug holds every row whose name contains SOMATROPIN/brand
tdrug['cls'] = (tdrug.drugname.fillna('') + ' ' + tdrug.prod_ai.fillna('')).map(drug_class)
def ids(cls, roles):
    if cls in ('DAILY_SOM', 'SERO_ZORB'):
        s = set(tdrug.loc[(tdrug.cls == cls) & tdrug.role_cod.isin(roles), 'primaryid'])
        s |= set(talld.loc[(talld.cls == cls) & talld.role_cod.isin(roles) & ~talld.primaryid.isin(tdrug.primaryid), 'primaryid'])
        return s & withpt
    s = set(talld.loc[(talld.cls == cls) & talld.role_cod.isin(roles), 'primaryid'])
    if cls == 'LONA': s |= set(lname.loc[lname.role_cod.isin(roles), 'primaryid'])
    return s & withpt
ALL = ['PS','SS','C','I']
SETS = {'L_PS': ids('LONA',['PS']), 'L_SUSP': ids('LONA',['PS','SS']), 'L_ANY': ids('LONA', ALL),
        'DS_PS': ids('DAILY_SOM',['PS']), 'SG_PS': ids('SOMATROGON',['PS']), 'SP_PS': ids('SOMAPACITAN',['PS']), 'SZ_PS': ids('SERO_ZORB',['PS']),
        'DS_ANY': ids('DAILY_SOM', ALL), 'SG_ANY': ids('SOMATROGON', ALL), 'SP_ANY': ids('SOMAPACITAN', ALL)}
OTHER_GH_ANY = SETS['DS_ANY'] | SETS['SG_ANY'] | SETS['SP_ANY'] | ids('SERO_ZORB', ALL)
R['counts'] = {k: len(v) for k, v in SETS.items()}
R['counts']['L_PS_with_other_GH_any_role'] = len(SETS['L_PS'] & OTHER_GH_ANY)
R['counts']['L_rows_by_name'] = (tdrug[tdrug.is_L].drugname + ' | ' + tdrug[tdrug.is_L].prod_ai).value_counts().to_dict()
R['counts']['L_reports_total_any_role_incl_no_pt'] = len(set(talld.loc[talld.cls == 'LONA', 'primaryid']) | set(lname.primaryid))
demo = tdemo.drop_duplicates('primaryid').set_index('primaryid')
AGEF = {'YR':1,'DEC':10,'MON':1/12,'WK':1/52,'DY':1/365,'HR':1/8760}
a_ = pd.to_numeric(demo.age, errors='coerce') * demo.age_cod.map(AGEF)
demo['age_y'] = a_.where((a_ >= 0) & (a_ <= 120))
demo['country'] = demo.occr_country.where(demo.occr_country.notna() & (demo.occr_country != 'None'), demo.reporter_country)
demo['country'] = demo.country.replace({'UNITED STATES':'US','None':np.nan})
demo['rcvd'] = demo.init_fda_dt.where(demo.init_fda_dt.notna() & (demo.init_fda_dt != 'None'), demo.fda_dt)
demo['year'] = pd.to_numeric(demo.rcvd.astype(str).str[:4], errors='coerce')
OCCP = {'MD':'Physician','HP':'Other health professional','PH':'Pharmacist','RN':'Nurse','CN':'Consumer','LW':'Lawyer'}
OUT = {'DE':'Death','LT':'Life-threatening','HO':'Hospitalization','DS':'Disability','CA':'Congenital anomaly','RI':'Required intervention','OT':'Other serious'}
AGEB = ([-0.1, 5.999, 11.999, 17.999, 64.999, 200], ['<6','6-11','12-17','18-64','>=65'])
def characterise(pids, label):
    d = demo.loc[sorted(pids)]; n = len(d); rows = []
    def add(var, cat, k): rows.append(dict(set=label, variable=var, category=cat, n=int(k), pct=round(100*k/n, 1) if n else None))
    for s, nm in [('F','Female'),('M','Male')]: add('Sex', nm, (d.sex == s).sum())
    add('Sex', 'Unknown/not reported', (~d.sex.isin(['F','M'])).sum())
    ag = pd.cut(d.age_y, AGEB[0], labels=AGEB[1])
    for c in AGEB[1]: add('Age (years)', c, (ag == c).sum())
    add('Age (years)', 'Unknown', ag.isna().sum())
    for c, nm in OCCP.items(): add('Reporter', nm, (d.occp_cod == c).sum())
    add('Reporter', 'Unknown/other', (~d.occp_cod.isin(OCCP)).sum())
    add('Country', 'United States', (d.country == 'US').sum()); add('Country', 'Other', (d.country.notna() & (d.country != 'US')).sum()); add('Country', 'Unknown', d.country.isna().sum())
    for c in ['EXP','PER','DIR']: add('Report type', {'EXP':'Expedited (15-day)','PER':'Periodic','DIR':'Direct'}[c], (d.rept_cod == c).sum())
    add('Report type', 'Other/unknown', (~d.rept_cod.isin(['EXP','PER','DIR'])).sum())
    for y in range(2021, 2027): add('Year of FDA receipt (retained version)', str(y), (d.quarter.str[:4] == str(y)).sum())
    o = toutc[toutc.primaryid.isin(pids)]
    for c, nm in OUT.items(): add('Outcome', nm, o.loc[o.outc_cod == c, 'primaryid'].nunique())
    add('Outcome', 'Any serious outcome', o.primaryid.nunique())
    st = dict(n=n, age_known=int(d.age_y.notna().sum()))
    if d.age_y.notna().any(): st.update(age_median=float(np.nanmedian(d.age_y)), age_q1=float(np.nanpercentile(d.age_y, 25)), age_q3=float(np.nanpercentile(d.age_y, 75)))
    return rows, st
rows = []
for k, lab in [('L_PS','Lonapegsomatropin (PS)'),('DS_PS','Daily somatropin (PS)'),('SG_PS','Somatrogon (PS)'),('SP_PS','Somapacitan (PS)')]:
    r_, R[f'{k}_age'] = characterise(SETS[k], lab); rows += r_
pd.DataFrame(rows).to_csv(f'{T}/table1_report_characteristics.csv', index=False)
R['L_PS_top_countries'] = demo.loc[sorted(SETS['L_PS'])].country.value_counts().head(8).to_dict()
R['L_PS_ped'] = int((demo.loc[sorted(SETS['L_PS'])].age_y < 18).sum()); R['L_PS_adult'] = int((demo.loc[sorted(SETS['L_PS'])].age_y >= 18).sum())

# ---------- quarterly pattern and 2026Q3 surge ----------
QS = sorted(pd.read_csv(f'{D}/quarter_stats.csv').quarter)
qt = pd.DataFrame({k: demo.loc[sorted(SETS[k])].quarter.value_counts().reindex(QS).fillna(0).astype(int) for k in ['L_PS','DS_PS','SG_PS','SP_PS']})
qs = pd.read_csv(f'{D}/quarter_stats.csv').set_index('quarter'); qt['N_all_reports'] = qs.N_all
nc_any_pids = set(treac.loc[treac.pt.map(lambda p: nonclinical(p)), 'primaryid'])
dl = demo.loc[sorted(SETS['L_PS'])]
qt['L_PS_with_nonclinical_PT'] = dl[dl.index.isin(nc_any_pids)].quarter.value_counts().reindex(QS).fillna(0).astype(int)
qt['L_PS_only_nonclinical_PTs'] = 0
clin_pids = set(treac.loc[~treac.pt.map(lambda p: nonclinical(p)), 'primaryid'])
qt['L_PS_only_nonclinical_PTs'] = dl[~dl.index.isin(clin_pids)].quarter.value_counts().reindex(QS).fillna(0).astype(int)
qt['L_PS_ped'] = dl[dl.age_y < 18].quarter.value_counts().reindex(QS).fillna(0).astype(int)
qt['L_PS_adult'] = dl[dl.age_y >= 18].quarter.value_counts().reindex(QS).fillna(0).astype(int)
qt['L_PS_age_unknown'] = dl[dl.age_y.isna()].quarter.value_counts().reindex(QS).fillna(0).astype(int)
qt['L_PS_consumer'] = dl[dl.occp_cod == 'CN'].quarter.value_counts().reindex(QS).fillna(0).astype(int)
qt['L_PS_US'] = dl[dl.country == 'US'].quarter.value_counts().reindex(QS).fillna(0).astype(int)
qt['L_PS_expedited'] = dl[dl.rept_cod == 'EXP'].quarter.value_counts().reindex(QS).fillna(0).astype(int)
qt['L_PS_periodic'] = dl[dl.rept_cod == 'PER'].quarter.value_counts().reindex(QS).fillna(0).astype(int)
qt['L_PS_direct'] = dl[dl.rept_cod == 'DIR'].quarter.value_counts().reindex(QS).fillna(0).astype(int)
qt['L_PS_initial_receipt_same_quarter'] = 0
idq = pd.to_datetime(dl.init_fda_dt.astype(str), format='%Y%m%d', errors='coerce').dt.to_period('Q').astype(str).str.replace('Q', 'Q')
qt['L_PS_initial_receipt_same_quarter'] = dl[idq.values == dl.quarter.values].quarter.value_counts().reindex(QS).fillna(0).astype(int)
qt.index.name = 'quarter'; qt.to_csv(f'{T}/table_quarterly_counts.csv')
prev4 = qt.loc[['2025Q3','2025Q4','2026Q1','2026Q2'], 'L_PS']
R['surge'] = dict(L_PS_2026Q3=int(qt.loc['2026Q3','L_PS']), prev4_mean=float(prev4.mean()), prev4=prev4.to_dict(),
                  ratio=float(qt.loc['2026Q3','L_PS']/prev4.mean()), max_before_2026Q3=int(qt.loc[QS[:-1],'L_PS'].max()),
                  quarter_of_max_before=str(qt.loc[QS[:-1],'L_PS'].idxmax()))
def comp_period(pids):
    d = demo.loc[sorted(pids)]; n = len(d)
    if n == 0: return {}
    has_nc = d.index.isin(nc_any_pids); only_nc = ~d.index.isin(clin_pids)
    return dict(n=n, ped=int((d.age_y < 18).sum()), adult=int((d.age_y >= 18).sum()), age_unknown=int(d.age_y.isna().sum()),
                consumer=int((d.occp_cod == 'CN').sum()), hcp=int(d.occp_cod.isin(['MD','HP','PH','RN']).sum()), us=int((d.country == 'US').sum()),
                expedited=int((d.rept_cod == 'EXP').sum()), periodic=int((d.rept_cod == 'PER').sum()), direct=int((d.rept_cod == 'DIR').sum()),
                with_nonclinical=int(has_nc.sum()), only_nonclinical=int(only_nc.sum()),
                serious=int(toutc[toutc.primaryid.isin(pids)].primaryid.nunique()),
                initial_fda_receipt_2026=int((pd.to_numeric(d.init_fda_dt.astype(str).str[:4], errors='coerce') == 2026).sum()),
                initial_fda_receipt_before_2025=int((pd.to_numeric(d.init_fda_dt.astype(str).str[:4], errors='coerce') < 2025).sum()))
q3 = set(dl.index[dl.quarter == '2026Q3']); pre = set(dl.index[dl.quarter != '2026Q3'])
R['surge']['composition_2026Q3'] = comp_period(q3); R['surge']['composition_before'] = comp_period(pre)
for nm, a_set, b_set in [('top_pts_2026Q3', q3, None), ('top_pts_before', pre, None)]:
    R['surge'][nm] = treac[treac.primaryid.isin(a_set)].pt.map(disp).value_counts().head(15).to_dict()
# PT shift 2026Q3 vs earlier (within lonapegsomatropin)
rq3 = treac[treac.primaryid.isin(q3)].pt.value_counts(); rpre = treac[treac.primaryid.isin(pre)].pt.value_counts()
shift = []
for pt in set(rq3.index) | set(rpre.index):
    a = int(rq3.get(pt, 0)); c = int(rpre.get(pt, 0))
    if a + c < 10: continue
    o = table2x2(a, len(q3)-a, c, len(pre)-c)
    shift.append(dict(pt_display=disp(pt), a_2026Q3=a, n_2026Q3=len(q3), c_before=c, n_before=len(pre), ROR=o['ROR'], lo=o['ROR_lo'], hi=o['ROR_hi'], nonclinical=nonclinical(pt)))
shift = pd.DataFrame(shift).sort_values('ROR', ascending=False); shift.to_csv(f'{T}/tableS_surge_pt_shift_2026Q3.csv', index=False)
R['surge']['pts_enriched_2026Q3'] = shift[(shift.lo > 1) & (shift.a_2026Q3 >= 10)].head(12)[['pt_display','a_2026Q3','c_before','ROR','lo','hi','nonclinical']].to_dict('records')
# initial receipt year of 2026Q3 reports (backlog vs new)
R['surge']['init_receipt_year_2026Q3'] = pd.to_numeric(dl.loc[sorted(q3)].init_fda_dt.astype(str).str[:4], errors='coerce').value_counts().sort_index().astype(int).to_dict()
R['surge']['caseversion_2026Q3_gt1'] = int((pd.to_numeric(dl.loc[sorted(q3)].caseversion, errors='coerce') > 1).sum())

# ---------- indications, dose, therapy duration ----------
lt = tdrug[tdrug.is_L & tdrug.primaryid.isin(SETS['L_ANY'])].copy()
ti = tindi.merge(lt[['primaryid','drug_seq']].drop_duplicates(), left_on=['primaryid','indi_drug_seq'], right_on=['primaryid','drug_seq'])
ti['indi_pt'] = ti.indi_pt.fillna('').str.upper()
R['indications_L_PS_top'] = ti[ti.primaryid.isin(SETS['L_PS'])].drop_duplicates(['primaryid','indi_pt']).indi_pt.map(disp).value_counts().head(15).to_dict()
R['routes_L_PS'] = lt[lt.primaryid.isin(SETS['L_PS'])].drop_duplicates('primaryid').route.fillna('Not reported').value_counts().head(6).to_dict()
th_all = tther.merge(lt[['primaryid','drug_seq']].drop_duplicates(), left_on=['primaryid','dsg_drug_seq'], right_on=['primaryid','drug_seq'])
th_all = th_all[th_all.primaryid.isin(SETS['L_PS'])]
dur_days = pd.to_numeric(th_all.dur, errors='coerce') * th_all.dur_cod.map({'YR':365,'MON':30.4,'WK':7,'DAY':1,'DY':1})
st_ = pd.to_datetime(th_all.start_dt.astype(str).where(th_all.start_dt.astype(str).str.len() == 8), format='%Y%m%d', errors='coerce')
en_ = pd.to_datetime(th_all.end_dt.astype(str).where(th_all.end_dt.astype(str).str.len() == 8), format='%Y%m%d', errors='coerce')
dur2 = (en_ - st_).dt.days
th_all['dur_d'] = dur_days.fillna(dur2)
long_ids = set(th_all.loc[th_all.dur_d >= 730, 'primaryid'])
R['long_term'] = dict(n_with_duration=int(th_all.dropna(subset=['dur_d']).primaryid.nunique()), n_ge_2y=len(long_ids),
                      top_pts_ge_2y=treac[treac.primaryid.isin(long_ids)].pt.map(disp).value_counts().head(10).to_dict())

# ---------- PT-level disproportionality ----------
def pt_table(pids, label, Nbg=N, bg_col='n_all'):
    r = treac[treac.primaryid.isin(pids)].drop_duplicates(['primaryid','pt'])
    cnt = r.pt.value_counts(); nT = len(pids); out = []
    for pt, a in cnt.items():
        nb = int(pc.loc[pt, bg_col]) if pt in pc.index else int(a)
        o = vs_all(int(a), nT, max(nb, int(a)), Nbg)
        out.append(dict(analysis=label, pt=pt, pt_display=disp(pt), a=int(a), nT=nT, pct_of_reports=round(100*a/nT, 1), nonclinical=nonclinical(pt),
                        nc_subcategory='; '.join(nc_subcats(pt)) if nonclinical(pt) else '', labeled=labeled(pt), indication_related=indication_related(pt),
                        groups='; '.join(g for g in GROUPS if in_group(pt, g)), **o))
    return pd.DataFrame(out).sort_values(['a','ROR'], ascending=[False, False])
pt_ps = pt_table(SETS['L_PS'], 'L_PS'); pt_ps.to_csv(f'{T}/table2_pt_disproportionality_L_PS_all.csv', index=False)
clin = pt_ps[~pt_ps.nonclinical]
sig = clin[clin.signal_primary].sort_values('ROR', ascending=False); sig.to_csv(f'{T}/table2_signals_L_PS_primary.csv', index=False)
ncs = pt_ps[pt_ps.nonclinical & pt_ps.signal_primary].sort_values('a', ascending=False); ncs.to_csv(f'{T}/tableS_nonclinical_signals_L_PS.csv', index=False)
R['pt_level'] = dict(n_distinct_pts=int(len(pt_ps)), n_clinical_pts=int(len(clin)), n_nonclinical_pts=int(pt_ps.nonclinical.sum()), n_clinical_pts_a_ge3=int((clin.a >= 3).sum()),
                     n_signal_primary=int(len(sig)), n_signal_highspec=int(sig.signal_highspec.sum()),
                     n_signal_primary_labeled=int(sig.labeled.sum()), n_signal_primary_unlabeled=int((~sig.labeled).sum()),
                     n_signal_primary_indication_related=int(sig.indication_related.sum()),
                     n_signal_primary_unlabeled_not_indication=int(((~sig.labeled) & (~sig.indication_related)).sum()),
                     n_sig_PRR=int(clin.sig_PRR.sum()), n_sig_ROR=int(clin.sig_ROR.sum()), n_sig_IC=int((clin.sig_IC & (clin.a >= 3)).sum()), n_sig_EB05=int(clin.sig_EB05.sum()),
                     n_nonclinical_signal_primary=int(len(ncs)),
                     top_reported_clinical=clin.head(20)[['pt_display','a','pct_of_reports']].to_dict('records'),
                     top_reported_nonclinical=pt_ps[pt_ps.nonclinical].head(15)[['pt_display','a','pct_of_reports']].to_dict('records'))
for lab in ['L_SUSP','L_ANY']:
    t = pt_table(SETS[lab], lab); t.to_csv(f'{T}/tableS_pt_disproportionality_{lab}.csv', index=False)
    c = t[~t.nonclinical]; R[f'pt_level_{lab}'] = dict(n_signal_primary=int(c.signal_primary.sum()), n_signal_highspec=int(c.signal_highspec.sum()))

# ---------- prespecified groups (incl. non-clinical subcategories) ----------
GLIST = list(GROUPS) + ['NC: any non-clinical'] + ['NC: '+n for n, _ in NC_SUB]
def ptin(p, g):
    if g.startswith('NC: '):
        if not nonclinical(p): return False
        return True if g == 'NC: any non-clinical' else g[4:] in nc_subcats(p)
    return in_group(p, g)
_memo = {}
def grp_ids(pids, g):
    r = treac[treac.primaryid.isin(pids)]
    key = g
    if key not in _memo: _memo[key] = {p for p in treac.pt.unique() if ptin(p, g)}
    return set(r.loc[r.pt.isin(_memo[key]), 'primaryid'])
def group_rows(pids, label, Nbg=N, col='n_all'):
    nT = len(pids); out = []
    for g in GLIST:
        a = len(grp_ids(pids, g)); bg = int(gb.loc[g, col]) if g in gb.index else a
        o = vs_all(a, nT, max(bg, a), Nbg)
        out.append(dict(analysis=label, group=g, a=a, nT=nT, pct=round(100*a/nT, 1) if nT else None, **o))
    return out
grp = pd.DataFrame(group_rows(SETS['L_PS'], 'L_PS')); grp.to_csv(f'{T}/table3_prespecified_groups_L_PS.csv', index=False)
R['groups_L_PS'] = grp.set_index('group')[['a','pct','ROR','ROR_lo','ROR_hi','IC025','EB05','signal_primary','signal_highspec']].to_dict('index')
comp_pts = {}
for g in GLIST:
    r = treac[treac.primaryid.isin(SETS['L_PS'])]
    comp_pts[g] = r[r.pt.isin(_memo.get(g, set()) or {p for p in treac.pt.unique() if ptin(p, g)})].drop_duplicates(['primaryid','pt']).pt.map(disp).value_counts().head(10).to_dict()
R['group_component_pts_L_PS'] = comp_pts

# ---------- active comparators ----------
KEY_PTS = ['Injection site pain','Injection site erythema','Injection site swelling','Injection site bruising','Injection site haemorrhage','Injection site pruritus',
           'Injection site reaction','Injection site atrophy','Lipoatrophy','Headache','Insulin-like growth factor increased','Insulin-like growth factor decreased',
           'Hyperglycaemia','Blood glucose increased','Scoliosis','Epiphysiolysis','Arthralgia','Pain in extremity','Oedema peripheral','Hypothyroidism',
           'Pyrexia','Vomiting','Urticaria','Rash','Hypersensitivity','Fatigue','Growth retardation','Dose omission issue','Device issue','Product dose omission issue',
           'Device malfunction','Product storage error','Incorrect dose administered','Underdose','Drug ineffective']
def comp_rows(A_, B_, label, bname):
    out = []
    ra = treac[treac.primaryid.isin(A_)]; rb = treac[treac.primaryid.isin(B_)]
    def one(name, a, c):
        o = table2x2(a, len(A_)-a, c, len(B_)-c)
        out.append(dict(analysis=label, item=name, a_L=a, n_L=len(A_), c_comp=c, n_comp=len(B_), comparator=bname,
                        ROR=o['ROR'], ROR_lo=o['ROR_lo'], ROR_hi=o['ROR_hi'], haldane=min(a, c, len(A_)-a, len(B_)-c) == 0,
                        higher=bool(a >= 3 and o['ROR_lo'] > 1), lower=bool(o['ROR_hi'] < 1)))
    for g in GLIST:
        S = _memo.get(g) or {p for p in treac.pt.unique() if ptin(p, g)}; _memo[g] = S
        one(g, ra.loc[ra.pt.isin(S), 'primaryid'].nunique(), rb.loc[rb.pt.isin(S), 'primaryid'].nunique())
    for pt in KEY_PTS:
        one('PT: '+pt, ra.loc[ra.pt == pt.upper(), 'primaryid'].nunique(), rb.loc[rb.pt == pt.upper(), 'primaryid'].nunique())
    return out
comp = pd.DataFrame(comp_rows(SETS['L_PS'], SETS['DS_PS'], 'L_PS vs DS_PS', 'Daily somatropin (GHD products)') +
                    comp_rows(SETS['L_PS'], SETS['SG_PS'], 'L_PS vs SG_PS', 'Somatrogon'))
comp.to_csv(f'{T}/table4_active_comparators.csv', index=False)
R['comparators'] = {lab: comp[comp.analysis == lab].set_index('item')[['a_L','n_L','c_comp','n_comp','ROR','ROR_lo','ROR_hi','higher','lower']].to_dict('index') for lab in comp.analysis.unique()}
for k in ['DS_PS','SG_PS']:
    t = pt_table(SETS[k], k); t.to_csv(f'{T}/tableS_pt_disproportionality_{k}.csv', index=False)
    R[f'{k}_signal_primary_n'] = int(t[~t.nonclinical].signal_primary.sum())
    R[f'{k}_groups'] = pd.DataFrame(group_rows(SETS[k], k)).set_index('group')[['a','pct','ROR','ROR_lo','ROR_hi','signal_primary']].round(3).to_dict('index')
def wilson(k, n):
    if n == 0: return (None, None)
    z = 1.96; p = k/n; den = 1 + z*z/n; c = (p + z*z/(2*n))/den; h = z*np.sqrt(p*(1-p)/n + z*z/(4*n*n))/den
    return (round(100*(c-h), 1), round(100*(c+h), 1))
sp_rows = []
for g in GLIST:
    S = _memo[g]
    for k, lab in [('L_PS','Lonapegsomatropin'),('DS_PS','Daily somatropin'),('SG_PS','Somatrogon'),('SP_PS','Somapacitan')]:
        n = len(SETS[k]); kk = treac[treac.primaryid.isin(SETS[k]) & treac.pt.isin(S)].primaryid.nunique(); lo, hi = wilson(kk, n)
        sp_rows.append(dict(group=g, product=lab, k=kk, n=n, pct=round(100*kk/n, 1) if n else None, lo=lo, hi=hi))
spd = pd.DataFrame(sp_rows); spd.to_csv(f'{T}/table4b_group_proportions_all_GH_products.csv', index=False)
R['somapacitan_n_PS'] = len(SETS['SP_PS']); R['somapacitan_ROR_computed'] = len(SETS['SP_PS']) >= 100
R['somapacitan_top_pts'] = treac[treac.primaryid.isin(SETS['SP_PS'])].pt.map(disp).value_counts().head(10).to_dict()

# ---------- pediatric vs adult ----------
ped = set(demo.index[demo.age_y < 18]); adu = set(demo.index[demo.age_y >= 18])
LP, LA = SETS['L_PS'] & ped, SETS['L_PS'] & adu
gp = pd.DataFrame(group_rows(LP, 'L_PS pediatric (<18) vs pediatric background', tot['N_ped'], 'n_ped'))
ga = pd.DataFrame(group_rows(LA, 'L_PS adult (>=18) vs adult background', tot['N_adult'], 'n_adult'))
m = gp.merge(ga, on='group', suffixes=('_ped', '_adult'))
def rorr(r):
    if r.a_ped < 3 or r.a_adult < 3: return pd.Series(dict(RORR=np.nan, RORR_lo=np.nan, RORR_hi=np.nan, p=np.nan))
    se1 = (np.log(r.ROR_hi_ped) - np.log(r.ROR_lo_ped))/(2*1.96); se2 = (np.log(r.ROR_hi_adult) - np.log(r.ROR_lo_adult))/(2*1.96)
    d = np.log(r.ROR_ped) - np.log(r.ROR_adult); se = np.sqrt(se1**2 + se2**2)
    return pd.Series(dict(RORR=np.exp(d), RORR_lo=np.exp(d-1.96*se), RORR_hi=np.exp(d+1.96*se), p=2*(1-norm.cdf(abs(d/se)))))
m = pd.concat([m, m.apply(rorr, axis=1)], axis=1)
m.to_csv(f'{T}/table5_pediatric_vs_adult_groups.csv', index=False)
R['age_strata'] = dict(n_ped=len(LP), n_adult=len(LA), n_unknown=len(SETS['L_PS']) - len(LP) - len(LA),
                       groups=m.set_index('group')[['a_ped','pct_ped','ROR_ped','ROR_lo_ped','ROR_hi_ped','signal_primary_ped','a_adult','pct_adult','ROR_adult','ROR_lo_adult','ROR_hi_adult','signal_primary_adult','RORR','RORR_lo','RORR_hi','p']].to_dict('index'))
for lab, S, Nb, col in [('ped', LP, tot['N_ped'], 'n_ped'), ('adult', LA, tot['N_adult'], 'n_adult')]:
    t = pt_table(S, f'L_PS_{lab}', Nb, col); t.to_csv(f'{T}/tableS_pt_disproportionality_L_PS_{lab}.csv', index=False)
    c = t[~t.nonclinical]; s = c[c.signal_primary].sort_values('ROR', ascending=False)
    R['age_strata'][f'{lab}_n_signal_primary'] = int(len(s)); R['age_strata'][f'{lab}_signals'] = s[['pt_display','a','ROR','ROR_lo','ROR_hi','IC025','EB05','labeled']].to_dict('records')
    R['age_strata'][f'{lab}_top_reported'] = t.head(10)[['pt_display','a','pct_of_reports','nonclinical']].to_dict('records')
dla = demo.loc[sorted(LA)]
R['adult_by_period'] = dict(before_2025Q3=int((dla.quarter < '2025Q3').sum()), from_2025Q3=int((dla.quarter >= '2025Q3').sum()),
                            by_quarter=dla.quarter.value_counts().sort_index().astype(int).to_dict())
dlp = demo.loc[sorted(LP)]; R['ped_age_bands'] = pd.cut(dlp.age_y, [-0.1, 5.999, 11.999, 17.999], labels=['<6','6-11','12-17']).value_counts().to_dict()

# ---------- sensitivity analyses ----------
SENS_T = {}
sens = []
def sens_summary(label, pids, Nbg=N, col='n_all'):
    t = pt_table(pids, label, Nbg, col) if len(pids) else pd.DataFrame()
    c = t[~t.nonclinical] if len(t) else t
    g = pd.DataFrame(group_rows(pids, label, Nbg, col)) if len(pids) else pd.DataFrame()
    row = dict(analysis=label, n_reports=len(pids), n_signal_primary=int(c.signal_primary.sum()) if len(c) else 0,
               n_signal_highspec=int(c.signal_highspec.sum()) if len(c) else 0)
    for _, x in g.iterrows():
        row[f"{x['group']}: a"] = x['a']; row[f"{x['group']}: ROR (95% CI)"] = f"{x['ROR']:.2f} ({x['ROR_lo']:.2f}-{x['ROR_hi']:.2f})"; row[f"{x['group']}: signal"] = x['signal_primary']
    sens.append(row); SENS_T[label] = t; return t
SENS_T = {}
sens_summary('Primary: PS', SETS['L_PS'])
sens_summary('S1: PS or SS', SETS['L_SUSP'])
sens_summary('S2: any role', SETS['L_ANY'])
us_ids = set(demo.index[demo.country == 'US']); hcp_ids = set(demo.index[demo.occp_cod.isin(['MD','HP','PH','RN'])])
sens_summary('S3: PS, US reports (US background)', SETS['L_PS'] & us_ids, tot['N_us'], 'n_us')
sens_summary('S4: PS, health-professional reporters (HCP background)', SETS['L_PS'] & hcp_ids, tot['N_hcp'], 'n_hcp')
dd = demo.loc[sorted(SETS['L_PS'])].copy()
ptset = treac[treac.primaryid.isin(SETS['L_PS'])].groupby('primaryid').pt.apply(lambda s: '|'.join(sorted(set(s))))
dd['ptset'] = ptset
dd['key'] = dd.sex.fillna('') + '_' + dd.age_y.round(0).astype(str) + '_' + dd.country.fillna('') + '_' + dd.event_dt.fillna('NA').astype(str) + '_' + dd.ptset.fillna('')
dupe = dd[dd.event_dt.notna() & (dd.event_dt.astype(str) != 'None') & dd.duplicated('key', keep='first')].index
R['S5_probable_duplicates_removed'] = int(len(dupe))
sens_summary('S5: PS, probable duplicates removed', SETS['L_PS'] - set(dupe))
sens_summary('S6: PS, excluding 2026Q3', SETS['L_PS'] - q3)
sens_summary('S7: PS pediatric vs pediatric background', LP, tot['N_ped'], 'n_ped')
pd.DataFrame(sens).to_csv(f'{T}/table6_sensitivity_analyses.csv', index=False)
rob = []
for _, x in sig.iterrows():
    r = dict(pt_display=x.pt_display, a=int(x.a), ROR=x.ROR, labeled=bool(x.labeled))
    for lab, tt_ in SENS_T.items():
        if lab.startswith('Primary') or lab.startswith('S7'): continue
        m_ = tt_[tt_.pt == x.pt] if len(tt_) else tt_
        r[lab.split(':')[0]] = bool(m_.signal_primary.iloc[0]) if len(m_) else False
    rob.append(r)
rob = pd.DataFrame(rob); rob.to_csv(f'{T}/table6b_signal_robustness.csv', index=False)
R['signal_robustness'] = rob.to_dict('records')
s8 = pd.DataFrame(comp_rows(SETS['L_PS'] - OTHER_GH_ANY, SETS['DS_PS'] - SETS['L_ANY'] - SETS['SG_ANY'] - SETS['SP_ANY'], 'S8: L_PS vs DS_PS excluding multi-GH reports', 'Daily somatropin') +
                  comp_rows(SETS['L_PS'] - OTHER_GH_ANY, SETS['SG_PS'] - SETS['L_ANY'] - SETS['DS_ANY'] - SETS['SP_ANY'], 'S8: L_PS vs SG_PS excluding multi-GH reports', 'Somatrogon'))
s8.to_csv(f'{T}/tableS_S8_comparators_excluding_multiGH.csv', index=False)
R['S8_n'] = dict(L=len(SETS['L_PS'] - OTHER_GH_ANY), DS=len(SETS['DS_PS'] - SETS['L_ANY'] - SETS['SG_ANY'] - SETS['SP_ANY']), SG=len(SETS['SG_PS'] - SETS['L_ANY'] - SETS['DS_ANY'] - SETS['SP_ANY']))
s9 = pd.DataFrame(comp_rows(SETS['L_PS'], SETS['DS_PS'] | SETS['SZ_PS'], 'S9: L_PS vs daily somatropin incl. Serostim/Zorbtive', 'Daily somatropin incl. Serostim/Zorbtive'))
s9.to_csv(f'{T}/tableS_S9_comparator_incl_serostim.csv', index=False); R['S9_n_comp'] = len(SETS['DS_PS'] | SETS['SZ_PS'])
# S11 (post hoc, documented deviation): comparators restricted to reports with >=1 clinical PT (daily somatropin reports are dominated by device/product PTs)
CLINREP = set(treac.loc[~treac.pt.map(lambda p: nonclinical(p)), 'primaryid'])
s11 = pd.DataFrame(comp_rows(SETS['L_PS'] & CLINREP, SETS['DS_PS'] & CLINREP, 'S12 (post hoc): L_PS vs DS_PS, reports with >=1 clinical PT', 'Daily somatropin') +
                   comp_rows(SETS['L_PS'] & CLINREP, SETS['SG_PS'] & CLINREP, 'S12 (post hoc): L_PS vs SG_PS, reports with >=1 clinical PT', 'Somatrogon'))
s11 = s11[~s11['item'].str.startswith('NC:')]
s11.to_csv(f'{T}/tableS_S12_comparators_clinical_reports_only.csv', index=False)
R['S12_n'] = dict(L=len(SETS['L_PS'] & CLINREP), DS=len(SETS['DS_PS'] & CLINREP), SG=len(SETS['SG_PS'] & CLINREP))
R['S12'] = {lab: s11[s11.analysis == lab].set_index('item')[['a_L','n_L','c_comp','n_comp','ROR','ROR_lo','ROR_hi','higher','lower']].to_dict('index') for lab in s11.analysis.unique()}
R['share_reports_only_nonclinical'] = {k: round(100*len(SETS[k] - CLINREP)/len(SETS[k]), 1) for k in ['L_PS','DS_PS','SG_PS','SP_PS']}
# S10: whole-database background (2010Q4-2026Q3, TESA-FAERS PT counts), PT level only
TD = f'{HERE}/../../tesamorelin-faers-2026-11/data/processed'
pcw = pd.read_parquet(f'{TD}/pt_counts_all_reports.parquet').set_index('pt'); Nw = json.load(open(f'{TD}/database_totals.json'))['N_reports_with_pt']
from lonastats import table2x2 as t22
s10 = []
for _, x in sig.iterrows():
    nb = int(pcw.loc[x.pt, 'n_all']) if x.pt in pcw.index else int(x.a)
    o = t22(int(x.a), int(x.nT) - int(x.a), nb - int(x.a), Nw - int(x.nT) - (nb - int(x.a)))
    s10.append(dict(pt_display=x.pt_display, a=int(x.a), ROR_window=x.ROR, ROR_wholeDB=o['ROR'], lo=o['ROR_lo'], hi=o['ROR_hi'], IC025=o['IC025'], signal=bool(x.a >= 3 and o['ROR_lo'] > 1 and o['IC025'] > 0)))
s10 = pd.DataFrame(s10); s10.to_csv(f'{T}/tableS_S10_whole_database_background.csv', index=False)
R['S10'] = dict(N_whole=int(Nw), n_primary_signals_retained=int(s10.signal.sum()) if len(s10) else 0, n_primary_signals=int(len(s10)))

# ---------- time to onset ----------
th = tther.merge(lt[['primaryid','drug_seq']].drop_duplicates(), left_on=['primaryid','dsg_drug_seq'], right_on=['primaryid','drug_seq'])
th = th[th.primaryid.isin(SETS['L_PS'])]
th['start'] = pd.to_datetime(th.start_dt.astype(str).where(th.start_dt.astype(str).str.len() == 8), format='%Y%m%d', errors='coerce')
evs = demo.event_dt.astype(str)
ev = pd.to_datetime(evs.where(evs.str.len() == 8), format='%Y%m%d', errors='coerce')
th['event'] = th.primaryid.map(ev)
th = th.dropna(subset=['start','event']).sort_values('start').drop_duplicates('primaryid', keep='first')
th['tto'] = (th.event - th.start).dt.days
n_dates = len(th); n_neg = int((th.tto < 0).sum()); n_long = int((th.tto > 5*365).sum()); th = th[(th.tto >= 0) & (th.tto <= 5*365)]
x = th.tto.values
def wb(xx, seed=1):
    if len(xx) < 10 or len(np.unique(xx)) < 3: return None
    c, loc, sc = weibull_min.fit(xx + 0.5, floc=0)
    rng = np.random.default_rng(seed); bs = []
    for _ in range(300):
        s = rng.choice(xx, len(xx), replace=True); bs.append(weibull_min.fit(s + 0.5, floc=0)[0])
    lo, hi = float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))
    return dict(shape=float(c), shape_lo=lo, shape_hi=hi if hi < 50 else None, scale=float(sc))
tto = dict(n_with_complete_dates=int(n_dates), n_negative_excluded=n_neg, n_over5y_excluded=n_long, n_valid=int(len(x)),
           median=float(np.median(x)) if len(x) else None, q1=float(np.percentile(x, 25)) if len(x) else None, q3=float(np.percentile(x, 75)) if len(x) else None,
           within_30_days=int((x <= 30).sum()), within_30_days_pct=round(100*(x <= 30).mean(), 1) if len(x) else None, weibull=wb(x))
by = {}
for g in ['Injection-site reactions','Hypersensitivity','Hyperglycaemia/type 2 diabetes','SCFE/osteonecrosis','Scoliosis','Intracranial hypertension','Headache','IGF-1 increased']:
    gpids = grp_ids(set(th.primaryid), g); xx = th.loc[th.primaryid.isin(gpids), 'tto'].values
    by[g] = dict(n=int(len(xx)), median=float(np.median(xx)) if len(xx) else None, q1=float(np.percentile(xx, 25)) if len(xx) else None,
                 q3=float(np.percentile(xx, 75)) if len(xx) else None, weibull=wb(xx) if len(xx) >= 10 else None)
tto['by_group'] = by; R['tto'] = tto
th[['primaryid','tto']].to_csv(f'{T}/tto_values_L_PS.csv', index=False)

# ---------- serious outcomes, deaths ----------
R['L_PS_deaths_n'] = int(toutc.loc[(toutc.outc_cod == 'DE') & toutc.primaryid.isin(SETS['L_PS']), 'primaryid'].nunique())
dpids = toutc.loc[(toutc.outc_cod == 'DE') & toutc.primaryid.isin(SETS['L_PS']), 'primaryid'].unique()
R['L_PS_death_pts'] = treac[treac.primaryid.isin(dpids)].pt.map(disp).value_counts().head(15).to_dict()
R['L_PS_death_ages'] = demo.loc[sorted(dpids)].age_y.tolist() if len(dpids) else []
# case-level descriptions for rare key groups (aggregate fields only)
cases = []
for g in ['Malignancy','CNS/pituitary tumour','SCFE/osteonecrosis','Intracranial hypertension','Pancreatitis','Hyperglycaemia/type 2 diabetes','Scoliosis']:
    for pid in sorted(grp_ids(SETS['L_PS'], g)):
        d = demo.loc[pid]; pts = treac[treac.primaryid == pid].pt.map(disp).tolist()
        cases.append(dict(group=g, age=None if pd.isna(d.age_y) else round(float(d.age_y), 1), sex=d.sex, quarter=d.quarter, reporter=d.occp_cod, country=d.country,
                          serious=bool(pid in set(toutc.primaryid)), group_pts='; '.join(p for p in pts if ptin(p.upper(), g)), n_pts=len(pts)))
pd.DataFrame(cases).to_csv(f'{T}/tableS_case_listing_key_groups_L_PS.csv', index=False)
R['key_group_case_summary'] = pd.DataFrame(cases).groupby('group').agg(n=('age', 'size'), n_age_lt18=('age', lambda s: int((s < 18).sum())), median_age=('age', 'median')).to_dict('index') if cases else {}
json.dump(R, open(f'{HERE}/../outputs/results.json', 'w'), indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o))
print(json.dumps({k: R[k] for k in ['counts','pt_level','surge','age_strata']}, indent=1, default=str)[:6000])
