"""LONAPEG-FAERS: derived numbers for the manuscript (no FAERS reprocessing). Rebuilds the L_PS case set exactly as 07_analysis.py
(asserted against results.json), then summarises the behavioural PT cluster at aggregate level, and computes percentages used in the text.
Writes outputs/derived_numbers.json. Aggregate values only (no case-level rows written)."""
import pandas as pd, numpy as np, json, os
from lonaconfig import LONA, drug_class
HERE = os.path.dirname(os.path.abspath(__file__)); D = f'{HERE}/../data/processed'; O = f'{HERE}/../outputs'; T = f'{O}/tables'
R = json.load(open(f'{O}/results.json')); V = json.load(open(f'{O}/review_results.json'))
tdrug = pd.read_parquet(f'{D}/target_tdrug.parquet'); tdemo = pd.read_parquet(f'{D}/target_tdemo.parquet'); treac = pd.read_parquet(f'{D}/target_treac.parquet')
talld = pd.read_parquet(f'{D}/target_talldrug.parquet'); toutc = pd.read_parquet(f'{D}/target_toutc.parquet')
for df in (tdrug, tdemo, treac, talld, toutc): df['primaryid'] = pd.to_numeric(df.primaryid).astype('int64')
withpt = set(treac.primaryid)
talld['cls'] = talld.dkey.map(drug_class)
isL = tdrug.drugname.str.contains(LONA) | tdrug.prod_ai.str.contains(LONA)
LPS = (set(talld.loc[(talld.cls == 'LONA') & (talld.role_cod == 'PS'), 'primaryid']) | set(tdrug.loc[isL & (tdrug.role_cod == 'PS'), 'primaryid'])) & withpt
assert len(LPS) == R['counts']['L_PS'], (len(LPS), R['counts']['L_PS'])
DISP = pd.read_csv(f'{D}/pt_display_names.csv').set_index('PT').display.to_dict()
demo = tdemo.drop_duplicates('primaryid').set_index('primaryid')
AGEF = {'YR':1,'DEC':10,'MON':1/12,'WK':1/52,'DY':1/365,'HR':1/8760}
a_ = pd.to_numeric(demo.age, errors='coerce') * demo.age_cod.map(AGEF); demo['age_y'] = a_.where((a_ >= 0) & (a_ <= 120))
rl = treac[treac.primaryid.isin(LPS)].copy(); rl['ptu'] = rl.pt.str.upper().str.strip()
BEH = ['AFFECTIVE DISORDER','AGGRESSION','MOOD SWINGS','CRYING','ABNORMAL BEHAVIOUR']
BEH_WIDE = BEH + ['IRRITABILITY','ANGER','MOOD ALTERED','EMOTIONAL DISORDER','PERSONALITY CHANGE','BEHAVIOUR DISORDER','ANXIETY','DEPRESSION','DEPRESSED MOOD','AGITATION','TEARFULNESS','EMOTIONAL DISTRESS','IMPULSIVE BEHAVIOUR','TEMPER TANTRUM','HOSTILITY']
def summ(ids):
    d = demo.loc[sorted(ids)]; ages = d.age_y.dropna()
    qs = d['quarter'].value_counts().to_dict() if 'quarter' in d else {}
    return dict(n_reports=int(len(ids)), n_age_known=int(len(ages)), n_age_lt18=int((ages < 18).sum()), age_median=float(ages.median()) if len(ages) else None,
                age_min=float(ages.min()) if len(ages) else None, age_max=float(ages.max()) if len(ages) else None,
                sex=d.sex.fillna('NR').value_counts().to_dict(), reporter=d.occp_cod.fillna('NR').value_counts().to_dict(), rept_cod=d.rept_cod.fillna('NR').value_counts().to_dict(),
                by_quarter={str(k): int(v) for k, v in sorted(qs.items())}, serious=int(len(set(ids) & set(toutc.primaryid))))
beh_ids = set(rl.loc[rl.ptu.isin(BEH), 'primaryid'])
out = {}
out['behavioural_cluster_5pt'] = summ(beh_ids)
out['behavioural_cluster_5pt']['pts_per_report'] = rl[rl.primaryid.isin(beh_ids) & rl.ptu.isin(BEH)].groupby('primaryid').ptu.nunique().value_counts().sort_index().to_dict()
out['behavioural_cluster_5pt']['co_reported_pts_top'] = rl[rl.primaryid.isin(beh_ids) & ~rl.ptu.isin(BEH)].ptu.map(lambda p: DISP.get(p, p.capitalize())).value_counts().head(12).to_dict()
present = sorted(set(rl.ptu) & set(BEH_WIDE))
out['behavioural_wide_pts_present'] = {DISP.get(p, p.capitalize()): int(rl[rl.ptu == p].primaryid.nunique()) for p in present}
out['behavioural_wide_n_reports'] = int(rl[rl.ptu.isin(BEH_WIDE)].primaryid.nunique())
gp = set(rl.loc[rl.ptu == 'GROWING PAINS', 'primaryid']); out['growing_pains'] = summ(gp)
ped = R['L_PS_ped']; ad = R['L_PS_adult']; n = R['counts']['L_PS']
t1 = pd.read_csv(f'{T}/table1_report_characteristics.csv'); L = t1[t1.set == 'Lonapegsomatropin (PS)']
g = lambda v, c: int(L[(L.variable == v) & (L.category == c)].n.iloc[0])
male, female = g('Sex', 'Male'), g('Sex', 'Female')
s = R['surge']; c3, cb = s['composition_2026Q3'], s['composition_before']
pct = lambda a, b: round(100 * a / b, 1)
out['pct'] = dict(age_unknown=pct(n - ped - ad, n), age_known=pct(ped + ad, n), male_of_known_sex=pct(male, male + female), us=pct(R['L_PS_top_countries']['US'], n),
    hcp=pct(g('Reporter', 'Physician') + g('Reporter', 'Other health professional') + g('Reporter', 'Pharmacist'), n),
    q3_share_of_all=pct(s['L_PS_2026Q3'], n), q3_periodic=pct(c3['periodic'], c3['n']), q3_consumer=pct(c3['consumer'], c3['n']), q3_only_nc=pct(c3['only_nonclinical'], c3['n']),
    q3_ped=pct(c3['ped'], c3['n']), before_only_nc=pct(cb['only_nonclinical'], cb['n']), before_periodic=pct(cb['periodic'], cb['n']), before_consumer=pct(cb['consumer'], cb['n']),
    before_age_unknown=pct(cb['age_unknown'], cb['n']), q3_age_unknown=pct(c3['age_unknown'], c3['n']), q3_init_2026=pct(c3['initial_fda_receipt_2026'], c3['n']),
    ped_of_age_known=pct(ped, ped + ad), any_nonclinical=R['groups_L_PS']['NC: any non-clinical']['pct'])
out['hcp_n'] = g('Reporter', 'Physician') + g('Reporter', 'Other health professional') + g('Reporter', 'Pharmacist')
out['S5_n_after'] = n - R['S5_probable_duplicates_removed']
out['signals_robust_all_S1_S6'] = [r['pt_display'] for r in R['signal_robustness'] if all(r[k] for k in ['S1','S2','S3','S4','S5','S6'])]
out['signals_lost_S6'] = [r['pt_display'] for r in R['signal_robustness'] if not r['S6']]
out['behav_signals_S4'] = {r['pt_display']: r['S4'] for r in R['signal_robustness'] if r['pt_display'] in ['Affective disorder','Aggression','Mood swings','Crying','Abnormal behaviour']}
out['behav_signals_S5'] = {r['pt_display']: r['S5'] for r in R['signal_robustness'] if r['pt_display'] in ['Affective disorder','Aggression','Mood swings','Crying','Abnormal behaviour']}
json.dump(out, open(f'{O}/derived_numbers.json', 'w'), indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o))
print(json.dumps(out, indent=1, default=str))
