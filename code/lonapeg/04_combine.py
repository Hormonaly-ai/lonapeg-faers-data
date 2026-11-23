"""LONAPEG-FAERS: combine per-quarter aggregates (window 2021Q3-2026Q3) into analysis files (data/processed)."""
import pandas as pd, glob, os, json
HERE=os.path.dirname(os.path.abspath(__file__)); OUT=f'{HERE}/../data/processed'; A=f'{OUT}/agg'
def cat(name): return pd.concat([pd.read_parquet(f) for f in sorted(glob.glob(f'{A}/*_{name}.parquet'))], ignore_index=True)
disp=cat('ptdisplay'); disp['PT']=disp.pt.str.upper()
disp=disp.groupby(['PT','pt']).k.sum().reset_index().sort_values('k',ascending=False).drop_duplicates('PT')[['PT','pt']].rename(columns={'pt':'display'})
disp.to_csv(f'{OUT}/pt_display_names.csv', index=False)
st=cat('stats'); st.to_csv(f'{OUT}/quarter_stats.csv', index=False)
cols=['n_all','n_us','n_hcp','n_ped','n_adult']
pc=cat('ptcounts').groupby('pt')[cols].sum().reset_index(); pc.to_parquet(f'{OUT}/pt_counts_all_reports.parquet', index=False)
gb=cat('gbg').groupby('group')[cols].sum().reset_index(); gb.to_csv(f'{OUT}/group_background_counts.csv', index=False)
marg=cat('marg').groupby('dkey').n_drug.sum().reset_index(); marg.to_parquet(f'{OUT}/mgps_drug_marginals.parquet', index=False)
pairs=pd.concat([pd.read_parquet(f) for f in sorted(glob.glob(f'{A}/*_pairs.parquet'))], ignore_index=True).groupby(['dkey','pt'],as_index=False).n.sum()
pairs.to_parquet(f'{OUT}/mgps_pairs_ps.parquet', index=False)
for n in ['tdrug','tdemo','treac','toutc','tther','tindi','talldrug']:
    cat(n).to_parquet(f'{OUT}/target_{n}.parquet', index=False)
s={('N_reports_with_pt' if k=='N_all' else k): int(st[k].sum()) for k in ['N_all','N_us','N_hcp','N_ped','N_adult','N_ps_reports']}
s.update(n_pairs=int(len(pairs)), n_drugs=int(len(marg)), n_pts=int(len(pc)), n_quarters=int(len(st)), quarters=f"{st.quarter.min()}-{st.quarter.max()}")
json.dump(s, open(f'{OUT}/database_totals.json','w'), indent=1); print(s)
