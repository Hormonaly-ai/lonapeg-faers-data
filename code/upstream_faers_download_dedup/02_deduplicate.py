"""Deduplication (protocol section 5): FAERS latest CASEVERSION per CASEID; legacy AERS latest ISR per CASE;
a case present in both systems keeps its FAERS version; CASEIDs in deleted-case files removed."""
import pandas as pd, glob, json, sys, os
P = sys.argv[1] if len(sys.argv)>1 else os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','data','processed','faers_reduced')
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', 'processed')
cols = ['primaryid','caseid','caseversion','fda_dt','quarter','system']
def load(f):
    x = pd.read_parquet(f, columns=cols)
    for c in ['primaryid','caseid','caseversion','fda_dt']: x[c] = pd.to_numeric(x[c], errors='coerce').astype('float64')
    x['quarter'] = x.quarter.astype('category'); x['system'] = x.system.astype('category')
    return x
parts = [load(f) for f in sorted(glob.glob(f'{P}/*_demo.parquet'))]
d = pd.concat(parts, ignore_index=True); del parts
d['quarter'] = d.quarter.astype('category'); d['system'] = d.system.astype('category')
n_rows = len(d); n_rows_sys = {k:int(v) for k,v in d.system.value_counts().items()}
d = d.dropna(subset=['primaryid','caseid'])
d['caseversion'] = d['caseversion'].fillna(0)
d = d.drop_duplicates('primaryid', keep='last')
n_unique_pid = len(d)
d['sysrank'] = (d.system == 'FAERS').astype(int); d['pid_abs'] = d.primaryid.abs()
d = d.sort_values(['caseid','sysrank','caseversion','fda_dt','pid_abs'])
latest = d.drop_duplicates('caseid', keep='last')
n_cases_both = int(d.groupby('caseid').system.nunique().gt(1).sum())
dels = pd.concat([pd.read_parquet(f) for f in glob.glob(f'{P}/*_deleted.parquet')])
dels = set(pd.to_numeric(dels['caseid'], errors='coerce').dropna().astype('int64'))
n_latest = len(latest)
latest = latest[~latest['caseid'].astype('int64').isin(dels)]
latest[['primaryid','caseid','caseversion','quarter','system']].astype({'primaryid':'int64','caseid':'int64','quarter':str,'system':str}).to_parquet(f'{OUT}/retained_reports.parquet', index=False)
summ = dict(demo_rows_all_quarters=int(n_rows), demo_rows_by_system=n_rows_sys, unique_primaryids=int(n_unique_pid), unique_caseids_latest_version=int(n_latest),
            caseids_in_both_systems=n_cases_both, deleted_caseids_listed=len(dels), removed_as_deleted=int(n_latest-len(latest)), retained_reports=int(len(latest)),
            retained_by_system={k:int(v) for k,v in latest.system.value_counts().items()}, n_quarters=int(d.quarter.nunique()))
json.dump(summ, open(f'{OUT}/dedup_summary.json','w'), indent=1); print(summ)
