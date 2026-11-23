"""LONAPEG-FAERS Part B: screening log, PRISMA counts, extraction tables (from full texts in review/fulltext; values transcribed by the reviewer
and traceable to the quoted passages in review/fulltext/_xml_safety_sections.md and the PDF text files)."""
import json, os, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); RV = f'{HERE}/../review'; O = f'{HERE}/../outputs'
S = json.load(open(f'{RV}/search/search_results.json'))
pm = {r['pmid'] for r in S['pubmed']}; pdoi = {r['doi'].lower() for r in S['pubmed'] if r['doi']}
extra = [r for r in S['epmc'] if (r['pmid'] not in pm) and (r['doi'].lower() not in pdoi or not r['doi'])]
n_pubmed, n_epmc = len(S['pubmed']), len(S['epmc']); n_dupes = n_epmc - len(extra); n_screened = n_pubmed + len(extra)
sought = [  # (key, id, design/stream, decision, reason)
 ('dka2026','PMID 42488375','Case report','Included',''),
 ('sat216','PMC12545226 (J Endocr Soc 2025;9 Suppl 1)','Case series (conference abstract; no fuller publication)','Included','abstract is the only publication'),
 ('jpem2025','PMID 39843388','Retrospective cohort (switch), n=14','Included',''),
 ('hyperglyc2023','PMID 37015214','Case report','Not retrieved','no UBC licence for Karger article (EZproxy, 23 Nov 2026); title/abstract only'),
 ('edm2026','PMID 41499205','Case report','Excluded','no adverse event described'),
 ('php1b2024','PMID 39193092','Case report','Excluded','no lonapegsomatropin safety data'),
 ('jei2026','PMID 42525213','Prospective cohort','Excluded','wrong intervention (somatrogon)'),
 ('insights2025','PMID 40702497','Registry interim analysis','Excluded','safety not reported by product'),
 ('sprogoe2018','PMID 29273485','Phase 1 RCT, healthy adults','Included',''),
 ('chatelain2017','PMID 28201598','Phase 2 RCT, children','Included',''),
 ('hoybye2017','PMID 28196799','Phase 2 RCT, adults','Included',''),
 ('thornton2021','PMID 34272849','Phase 3 RCT heiGHt','Included',''),
 ('maniatis2022f','PMID 35263755','Phase 3 single-arm fliGHt','Included',''),
 ('maniatis2022e','PMID 35428884','Phase 3 extension enliGHten (interim)','Included','companion report of enliGHten'),
 ('enlighten_final','PMID 40049149','Phase 3 extension enliGHten (final)','Included',''),
 ('foresight','PMID 41420532 (+ correction PMID 42842218)','Phase 3 RCT foresiGHt (adults)','Included',''),
 ('adult_ole','PMID 42855131','foresiGHt open-label extension','Included',''),
 ('china','PMID 41758754','Phase 3 RCT, China','Included',''),
 ('turner','PMID 42929032','Phase 2 RCT, Turner syndrome (New InsiGHTS)','Included',''),
 ('coach','PMID 42144862','Phase 2 single-arm, achondroplasia (with navepegritide)','Included',''),
 ('mygind2026','PMID 42768234','Pooled immunogenicity analysis (heiGHt, fliGHt, enliGHten)','Included','companion analysis'),
]
sd = pd.DataFrame(sought, columns=['key','identifier','design','decision','reason']); sd.to_csv(f'{RV}/extraction/fulltext_decisions.csv', index=False)
inc = sd[sd.decision == 'Included']
P = dict(records_pubmed=n_pubmed, records_epmc=n_epmc, duplicates_removed=n_dupes, records_screened=n_screened,
         records_excluded_title_abstract=n_screened - len(sd), reports_sought=len(sd), reports_not_retrieved=int((sd.decision == 'Not retrieved').sum()),
         reports_assessed=int((sd.decision != 'Not retrieved').sum()), reports_excluded=int((sd.decision == 'Excluded').sum()),
         excluded_reasons=sd[sd.decision == 'Excluded'].reason.value_counts().to_dict(), reports_included=len(inc),
         stream1_reports=3, stream2_reports=len(inc) - 3, stream2_unique_studies=11, regulatory_documents=7,
         ctgov_records=14, ctgov_new_studies=0, search_date='2026-11-23')
# Stream 2 extraction (lonapegsomatropin arms; numbers from full texts)
T2 = [
 dict(study='Phase 1 (healthy men)', ref='sprogoe2018', population='Healthy adult men', N_lona=28, duration='single dose', teae='reported by cohort; headache most common', sae='0', disc_ae='NR', deaths='0', isr='no difference vs Omnitrope or placebo; no nodules or lipoatrophy', igf='dose-dependent IGF-1 rise', glucose='NR', other='decreased neutrophil count across all groups incl. placebo', ada='NR'),
 dict(study='Phase 2 (children)', ref='chatelain2017', population='Prepubertal GHD, treatment-naive', N_lona=None, duration='26 wk', teae='29/53 (54.7%) incl. Genotropin arm', sae='1/53 (inguinal hernia, unrelated)', disc_ae='NR', deaths='0', isr='25 subjects, 141 events (pain most common); no nodules or lipoatrophy', igf='dose-dependent', glucose='stable', other='', ada='1 low-titre non-neutralising'),
 dict(study='Phase 2 (adults)', ref='hoybye2017', population='Adult GHD', N_lona=29, duration='4 weekly doses', teae='22/29 (75.9%) across 3 doses', sae='1 participant (adrenal crisis, pleuritic chest pain; unrelated)', disc_ae='3/29', deaths='0', isr='application-site erythema/induration 2/10 (lowest dose); no lipoatrophy or nodules', igf='dose-dependent', glucose='stable', other='affect lability, depressed mood and insomnia each 2/10 at lowest dose', ada='none'),
 dict(study='heiGHt (phase 3 RCT)', ref='thornton2021', population='Prepubertal GHD, treatment-naive', N_lona=105, duration='52 wk', teae='81/105 (77.1%) vs 39/56 (69.6%) daily', sae='1 (appendicitis, unrelated)', disc_ae='0', deaths='0', isr='lipoatrophy 1 (same-site injections), urticaria 1', igf='dose reduction for raised IGF-1 2/105 (1.9%)', glucose='HbA1c and fasting glucose stable', other='secondary hypothyroidism 6.7% vs 7.1%; adrenal insufficiency 1.9% vs 1.8%', ada='7/105 (6.7%) low-titre, non-neutralising'),
 dict(study='fliGHt (phase 3)', ref='maniatis2022f', population='Children switching from daily GH', N_lona=146, duration='26 wk', teae='83/146 (57%); related 4%', sae='2 in 1 participant (AV block, chest pain; unrelated)', disc_ae='0', deaths='0', isr='>=96.6% no redness/bruising/swelling/pain at 2 h', igf='NR here', glucose='HbA1c 5.2% stable', other='', ada='3% low-titre anti-hGH; no anti-PEG'),
 dict(study='enliGHten (extension, up to 6 y)', ref='enlighten_final', population='Children from heiGHt/fliGHt', N_lona=298, duration='up to 6 y (mean 4.1 y)', teae='226/298 (75.8%)', sae='35 SAEs in 21 (7%); all unrelated', disc_ae='0', deaths='0', isr='12 ISRs in 8 participants (related); local reactions 7.1% of vial/syringe vs 1.3% of auto-injector injections', igf='dose reductions in 89 (29.9%) (interim report, Maniatis 2022)', glucose='HbA1c stable', other='scoliosis 2 (related, mild/moderate); no intracranial hypertension or pancreatitis (SCFE not mentioned; corrected at audit 23 Nov 2026)', ada='25/298 (8.4%) any ADA; anti-mPEG 2/298 (0.7%) treatment-emergent; no neutralising'),
 dict(study='foresiGHt (phase 3 RCT, adults)', ref='foresight', population='Adult GHD', N_lona=89, duration='38 wk', teae='64/89 (71.9%) vs 65.5% placebo, 73.3% somatropin', sae='4 (4.5%); 1 related (hyponatraemia)', disc_ae='1 (1.1%)', deaths='0', isr='4 (4.5%) vs 4.8% placebo, 5.8% somatropin', igf='titrated', glucose='no new-onset diabetes with lonapegsomatropin', other='arthralgia 9.0%, headache 7.9%', ada='see extension'),
 dict(study='foresiGHt extension (adults)', ref='adult_ole', population='Adult GHD', N_lona=236, duration='38 + 52 wk', teae='168/236 (71.2%)', sae='20 (8.5%); all unrelated', disc_ae='9 (3.8%)', deaths='NR', isr='12 (5.1%)', igf='titrated', glucose='HbA1c reported over time', other='basal cell carcinoma 1 and progression of 2 pre-existing benign neoplasms (unrelated); no recurrence on cranial MRI in childhood cancer survivors', ada='anti-hGH 3.4%, anti-lonapegsomatropin 3.4%, anti-PEG 0.9%'),
 dict(study='China phase 3 RCT', ref='china', population='Prepubertal GHD, China', N_lona=100, duration='52 wk', teae='98/100 (98.0%) vs 50/53 (94.3%)', sae='4 (4.0%); 1 possibly related (subcutaneous abscess)', disc_ae='2 (allergic dermatitis with abscess; gynaecomastia)', deaths='0', isr='5 (5.0%) vs 4 (7.5%)', igf='NR here', glucose='HbA1c stable; no new diabetes', other='', ada='NR here'),
 dict(study='New InsiGHTS (Turner syndrome)', ref='turner', population='Prepubertal Turner syndrome', N_lona=37, duration='52 wk', teae='see Table 2 of source (ISR 15/37, 40.5%)', sae='0', disc_ae='0', deaths='0', isr='15/37 (40.5%) vs 5/12 (41.7%) somatropin (combined term)', igf='dose-dependent', glucose='NR here', other='primary hypothyroidism 1; no intracranial hypertension or SCFE', ada='NR here'),
 dict(study='COACH (achondroplasia, with navepegritide)', ref='coach', population='Achondroplasia, 2-11 y', N_lona=21, duration='52 wk', teae='treatment-related 33.3% (lonapegsomatropin)', sae='2 in 1 child (hydrocele, viral RTI; unrelated)', disc_ae='0', deaths='0', isr='transient, mild', igf='NR here', glucose='NR here', other='dosing interrupted 10 weeks in 1 child for macrocephaly work-up', ada='NR here'),
]
t2 = pd.DataFrame(T2); t2.to_csv(f'{RV}/extraction/stream2_trials.csv', index=False)
T1 = [
 dict(ref='dka2026', country='USA', age_sex='10 y 10 mo, male', indication='GHD after suprasellar germ cell tumour (chemotherapy, craniospinal irradiation); central hypothyroidism, adrenal insufficiency, AVP deficiency; obesity (BMI 35.3)',
      dose='0.23 mg/kg/week', event='Diabetic ketoacidosis (glucose 333 mg/dL, bicarbonate 6 mmol/L; autoantibodies negative; C-peptide 26.6 ng/mL; IGF-1 333 ng/mL)', latency='first dose; vomiting day 3, DKA day 4',
      cofactors='stress-dose hydrocortisone for vomiting; aripiprazole; obesity; maternal type 2 diabetes', dechallenge='GH stopped; normoglycaemic thereafter', rechallenge='No', outcome='Recovered', naranjo='4 (possible)', jbi='7/8 yes'),
 dict(ref='sat216', country='USA (single paediatric centre)', age_sex='7 patients, mean 8.2 y (4.9-13.8), 6 male', indication='GHD (3 treatment-naive, 4 switched from daily GH)', dose='NR', event='lipoatrophy 3; marked behavioural change 1; hypersensitivity 3 (oedema, rash, pruritus), 2 systemic requiring urgent care',
      latency='NR', cofactors='NR', dechallenge='mostly resolved after stopping or with site rotation', rechallenge='NR', outcome='2 switched to another weekly GH, 3 to daily GH, 1 continued, 1 lost to follow-up', naranjo='not assessable (aggregate abstract)', jbi='case series 4/10 yes (abstract only)'),
 dict(ref='jpem2025', country='USA (single centre)', age_sex='14 children, median 9 y', indication='GHD, switched from daily somatropin', dose='median 0.23 mg/kg/week', event='IGF-1 +2.3 SDS after switch; 5 needed immediate dose reduction (4 further); leg pain 4 (SCFE excluded); headache/migraine 2 (normal optic nerve exam); ISRs 2; 1 death from septic shock (drug held; not attributed)',
      latency='leg pain within several months; headaches up to 1 y', cofactors='NR', dechallenge='leg pain resolved without intervention', rechallenge='NR', outcome='no serious adverse effects attributed', naranjo='not applicable (cohort)', jbi='case series 8/10 yes'),
]
t1 = pd.DataFrame(T1); t1.to_csv(f'{RV}/extraction/stream1_cases.csv', index=False)
REG = ['FDA BLA 761177 approval package (2021)', 'FDA BLA 761177/S-008 approval letter (adult GHD, 25 Jul 2025)', 'Skytrofa US prescribing information (openFDA version 14, effective 21 Oct 2026)',
       'EMA EPAR public assessment report (EMA/706519/2021) incl. RMP safety concerns', 'EMA PSUSA/00010969/202202 scientific conclusions (anaphylaxis/angioedema added)',
       'EMA PSUSA/00010969/202402 scientific conclusions (osteonecrosis as potential risk)', 'EMA product information (last updated 25 Sep 2026)']
V = dict(prisma=P, stream1=T1, stream2=T2, regulatory=REG,
         stream2_total_lona_exposed_main_trials=int(sum(x['N_lona'] or 0 for x in T2 if x['ref'] in ['thornton2021','maniatis2022f','foresight','china','turner','coach','hoybye2017','sprogoe2018'])))
json.dump(V, open(f'{O}/review_results.json', 'w'), indent=1); print(json.dumps(P, indent=1))
