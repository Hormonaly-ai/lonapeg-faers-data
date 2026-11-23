"""Reduce one FAERS (2012Q4+) or legacy AERS (2010Q4-2012Q3) quarterly ASCII zip to parquet (TESA-FAERS).
Usage: python 01_process_quarter.py /path/{faers|aers}_ascii_YYYYQn.zip OUTDIR
Keeps for all reports: DEMO key fields, DRUG (primaryid, role_cod, dkey, target flags), REAC (primaryid, pt), OUTC,
primaryids with an HIV-related indication (any drug); and full DRUG/THER/INDI detail rows for target drugs
(tesamorelin and somatropin). Legacy AERS: ISR -> primaryid (stored negative to avoid collisions), CASE -> caseid,
GNDR_COD -> sex, no PROD_AI (empty)."""
import sys, zipfile, re, hashlib, json, os
import pandas as pd
zp, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)
m = re.search(r'(\d{4})[qQ](\d)', os.path.basename(zp)); tag = f"{m.group(1)}Q{m.group(2)}"
legacy = os.path.basename(zp).lower().startswith('aers')
z = zipfile.ZipFile(zp); names = z.namelist()
def find(prefix):
    c = [n for n in names if re.search(rf'(^|/){prefix}\d\dQ\d(_new)?\.txt$', n, re.I)]
    return c[0]
def clean(c): return c.replace('\ufeff','').replace('ï»¿','').strip().lower()
def read(prefix, usecols):
    with z.open(find(prefix)) as f:
        df = pd.read_csv(f, sep='$', dtype=str, usecols=lambda c: clean(c) in usecols,
                         encoding='latin-1', on_bad_lines='skip', quoting=3, low_memory=False, index_col=False)
    df.columns = [clean(c) for c in df.columns]
    return df
TES = re.compile(r'TESAMORELIN|EGRIFTA|\bTH[- ]?9507\b', re.I)
SOM = re.compile(r'SOMATROPIN|GENOTROPIN|HUMATROPE|NORDITROPIN|NUTROPIN|OMNITROPE|SAIZEN|SEROSTIM|ZOMACTON|ZORBTIVE|TEV-TROPIN|ACCRETROPIN|VALTROPIN', re.I)
if legacy:
    demo = read('DEMO', {'isr','case','i_f_cod','foll_seq','event_dt','fda_dt','rept_cod','age','age_cod','gndr_cod','wt','wt_cod','occp_cod','reporter_country','rept_dt','mfr_dt'})
    demo = demo.rename(columns={'isr':'primaryid','case':'caseid','gndr_cod':'sex','i_f_cod':'i_f_code'})
    demo['caseversion'] = None; demo['occr_country'] = None; demo['init_fda_dt'] = None
    drug = read('DRUG', {'isr','drug_seq','role_cod','drugname','route','dose_vbm','dechal','rechal','lot_num','nda_num'}).rename(columns={'isr':'primaryid'})
    drug['prod_ai'] = ''
    reac = read('REAC', {'isr','pt'}).rename(columns={'isr':'primaryid'})
    outc = read('OUTC', {'isr','outc_cod'}).rename(columns={'isr':'primaryid'})
    ther = read('THER', {'isr','drug_seq','start_dt','end_dt','dur','dur_cod'}).rename(columns={'isr':'primaryid','drug_seq':'dsg_drug_seq'})
    indi = read('INDI', {'isr','drug_seq','indi_pt'}).rename(columns={'isr':'primaryid','drug_seq':'indi_drug_seq'})
    for df in (demo, drug, reac, outc, ther, indi):
        df['primaryid'] = '-' + df['primaryid'].astype(str).str.strip()
else:
    demo = read('DEMO', {'primaryid','caseid','caseversion','i_f_code','i_f_cod','event_dt','init_fda_dt','fda_dt','rept_cod','age','age_cod','sex','gndr_cod','wt','wt_cod','occp_cod','reporter_country','occr_country','rept_dt','mfr_dt'})
    demo = demo.rename(columns={'gndr_cod':'sex','i_f_cod':'i_f_code'})
    drug = read('DRUG', {'primaryid','drug_seq','role_cod','drugname','prod_ai','route','dose_vbm','dose_amt','dose_unit','dose_freq','dechal','rechal','lot_num','nda_num'})
    if 'prod_ai' not in drug: drug['prod_ai'] = ''
    reac = read('REAC', {'primaryid','pt','drug_rec_act'})
    outc = read('OUTC', {'primaryid','outc_cod','outc_code'}).rename(columns={'outc_code':'outc_cod'})
    ther = read('THER', {'primaryid','dsg_drug_seq','start_dt','end_dt','dur','dur_cod'})
    indi = read('INDI', {'primaryid','indi_drug_seq','indi_pt'})
demo['quarter'] = tag; demo['system'] = 'AERS' if legacy else 'FAERS'
drug['drugname'] = drug['drugname'].fillna('').str.upper().str.strip()
drug['prod_ai'] = drug['prod_ai'].fillna('').str.upper().str.strip()
drug['role_cod'] = drug['role_cod'].fillna('').str.upper().str.strip()
reac['pt'] = reac['pt'].fillna('').str.strip()
drug['is_tes'] = drug['drugname'].str.contains(TES) | drug['prod_ai'].str.contains(TES)
drug['is_som'] = drug['drugname'].str.contains(SOM) | drug['prod_ai'].str.contains(SOM)
drug['is_serostim'] = drug['drugname'].str.contains('SEROSTIM')
tdrug = drug[drug.is_tes | drug.is_som].copy()
tpids = set(tdrug['primaryid'])
drug['dkey'] = drug['prod_ai'].where(drug['prod_ai'] != '', drug['drugname'])
drug[['primaryid','role_cod','dkey','is_tes','is_som','is_serostim']].to_parquet(f'{out}/{tag}_drug.parquet', index=False)
tdrug.to_parquet(f'{out}/{tag}_tdrug.parquet', index=False)
demo.to_parquet(f'{out}/{tag}_demo.parquet', index=False)
reac[['primaryid','pt']].to_parquet(f'{out}/{tag}_reac.parquet', index=False)
outc[['primaryid','outc_cod']].to_parquet(f'{out}/{tag}_outc.parquet', index=False)
ther[ther['primaryid'].isin(tpids)].to_parquet(f'{out}/{tag}_tther.parquet', index=False)
indi[indi['primaryid'].isin(tpids)].to_parquet(f'{out}/{tag}_tindi.parquet', index=False)
hiv = indi[indi['indi_pt'].fillna('').str.upper().str.contains(r'\bHIV\b|HUMAN IMMUNODEFICIENCY|ACQUIRED IMMUNODEFICIENCY|\bAIDS\b', regex=True)]
pd.DataFrame({'primaryid': sorted(set(hiv['primaryid']))}).to_parquet(f'{out}/{tag}_hivindi.parquet', index=False)
dels = []
for n in names:
    if re.search(r'delete', n, re.I) and n.lower().endswith('.txt'):
        txt = z.read(n).decode('latin-1')
        dels += [l.strip() for l in txt.splitlines() if l.strip().isdigit()]
pd.DataFrame({'caseid': dels}).to_parquet(f'{out}/{tag}_deleted.parquet', index=False)
md5 = hashlib.md5(open(zp,'rb').read()).hexdigest()
log = dict(quarter=tag, system='AERS' if legacy else 'FAERS', zip=os.path.basename(zp), zip_bytes=os.path.getsize(zp), md5=md5,
           n_demo=len(demo), n_drug=len(drug), n_reac=len(reac), n_tes_rows=int(drug.is_tes.sum()), n_som_rows=int(drug.is_som.sum()),
           n_hiv_indi_reports=int(hiv.primaryid.nunique()), n_deleted=len(dels))
with open(f'{out}/process_log.jsonl','a') as f: f.write(json.dumps(log)+'\n')
print(log)
