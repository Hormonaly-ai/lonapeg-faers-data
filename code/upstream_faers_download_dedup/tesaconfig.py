"""Prespecified definitions for TESA-FAERS (protocol v1.0, OSF 10.17605/OSF.IO/H5PBN)."""
import re
TES = re.compile(r'TESAMORELIN|EGRIFTA|\bTH[- ]?9507\b', re.I)
SOM = re.compile(r'SOMATROPIN|GENOTROPIN|HUMATROPE|NORDITROPIN|NUTROPIN|OMNITROPE|SAIZEN|SEROSTIM|ZOMACTON|ZORBTIVE|TEV-TROPIN|ACCRETROPIN|VALTROPIN', re.I)
SERO = re.compile(r'SEROSTIM', re.I)

def formulation(drugname):
    d = (drugname or '').upper()
    if 'EGRIFTA' in d:
        if re.search(r'\bWR\b', d): return 'Egrifta WR'
        if re.search(r'\bSV\b', d): return 'Egrifta SV'
        return 'Egrifta (unspecified/1 mg)'
    if re.search(r'COMPOUND|RESEARCH|PEPTIDE', d): return 'Unbranded tesamorelin (compounded/research flag)'
    return 'Unbranded tesamorelin'

def _exact(*xs): return {x.upper() for x in xs}
GROUPS = {
 'Hypersensitivity': dict(exact=_exact('Hypersensitivity','Drug hypersensitivity','Anaphylactic reaction','Anaphylactic shock','Anaphylactoid reaction',
    'Urticaria','Rash','Rash pruritic','Rash erythematous','Rash generalised','Rash macular','Rash papular','Pruritus','Pruritus generalised',
    'Erythema','Angioedema','Swelling face','Lip swelling','Flushing')),
 'Injection-site reactions': dict(startswith=['INJECTION SITE']),
 'Fluid retention/musculoskeletal': dict(exact=_exact('Oedema peripheral','Oedema','Peripheral swelling','Fluid retention','Swelling','Generalised oedema',
    'Arthralgia','Pain in extremity','Myalgia','Joint swelling','Musculoskeletal stiffness','Musculoskeletal pain','Paraesthesia','Hypoaesthesia')),
 'Carpal tunnel syndrome': dict(exact=_exact('Carpal tunnel syndrome')),
 'Hyperglycaemia/diabetes': dict(exact=_exact('Hyperglycaemia','Blood glucose increased','Diabetes mellitus','Type 2 diabetes mellitus',
    'Diabetes mellitus inadequate control','Glycosylated haemoglobin increased','Glucose tolerance impaired','Impaired fasting glucose','Insulin resistance','Diabetic ketoacidosis')),
 'IGF-1 increased': dict(exact=_exact('Insulin-like growth factor increased','Blood growth hormone increased')),
 'Malignancy/neoplasm': dict(contains=['MALIGNAN','CARCINOMA','CANCER','LYMPHOMA','LEUKAEMIA','SARCOMA','MELANOMA','METASTA','MYELOMA'], neoplasm_not_benign=True),
 'Intracranial hypertension': dict(exact=_exact('Idiopathic intracranial hypertension','Intracranial pressure increased','Benign intracranial hypertension','Papilloedema')),
 'Cardiovascular': dict(exact=_exact('Myocardial infarction','Acute myocardial infarction','Cerebrovascular accident','Ischaemic stroke','Transient ischaemic attack',
    'Cardiac failure','Cardiac failure congestive','Atrial fibrillation','Arrhythmia','Tachycardia','Palpitations','Chest pain','Hypertension',
    'Blood pressure increased','Deep vein thrombosis','Pulmonary embolism')),
 'Retinopathy': dict(exact=_exact('Diabetic retinopathy','Retinopathy')),
 'Pregnancy exposure': dict(exact=_exact('Exposure during pregnancy','Maternal exposure during pregnancy','Pregnancy','Abortion spontaneous')),
}
def in_group(pt, g):
    p = pt.upper(); d = GROUPS[g]
    if p in d.get('exact', set()): return True
    if any(s in p for s in d.get('contains', [])): return True
    if d.get('neoplasm_not_benign') and 'NEOPLASM' in p and 'BENIGN' not in p: return True
    return any(p.startswith(s) for s in d.get('startswith', []))

# Non-clinical PTs excluded from signal ranking (protocol section 7); reported descriptively (product/dose issues)
NONCLIN = re.compile(r'OFF LABEL|PRODUCT |^PRODUCT\b|PRODUCT$|DRUG INEFFECTIVE|WRONG TECHNIQUE|NO ADVERSE EVENT|INTENTIONAL PRODUCT|DEVICE|DOSE OMISSION|'
                     r'INAPPROPRIATE SCHEDULE|INCORRECT DOSE|EXTRA DOSE|WRONG (DRUG|PATIENT|PRODUCT|DOSE|ROUTE)|MEDICATION ERROR|EXPIRED PRODUCT|'
                     r'THERAPEUTIC RESPONSE (UNEXPECTED|DECREASED)|DRUG EFFECT (LESS|DECREASED|DELAYED|INCOMPLETE)|TREATMENT NONCOMPLIANCE|NEEDLE ISSUE|'
                     r'UNDERDOSE|ACCIDENTAL EXPOSURE|CIRCUMSTANCE OR INFORMATION|ADVERSE EVENT$|^ILLNESS$|INSURANCE|PRESCRIBED|DRUG DOSE OMISSION|'
                     r'INTERCEPTED|SUSPECTED COUNTERFEIT|DRUG DIVERSION|NO REACTION|UNEVALUABLE EVENT|CONDITION AGGRAVATED|DRUG ABUSE|'
                     r'SYRINGE ISSUE|INJURY ASSOCIATED WITH DEVICE|LIQUID PRODUCT|DRUG ADMINISTERED|LACK OF(?! SPONTANEOUS SPEECH)|DRUG DISPENSING|UNINTENTIONAL|'
                     r'RECONSTITUTION|POOR QUALITY (?!SLEEP)|DRUG DOSE TITRATION|DOSE CALCULATION|DRUG PRESCRIBING|INCORRECT ROUTE|'
                     r'INCORRECT PRODUCT|MISSED DOSE|INCORRECT (DRUG )?DOSAGE|ADMINISTRATION ERROR|TITRATION ERROR|DRUG DELIVERY SYSTEM|ADVERSE DRUG REACTION$|(?<!SEVERE CUTANEOUS )ADVERSE REACTION$|THERAPY INTERRUPTED|THERAPY CESSATION|'
                     r'THERAPY NON-RESPONDER|TREATMENT FAILURE|DRUG INTOLERANCE|INTENTIONAL DOSE|OVERDOSE|\bVIAL\b|DRUG STORAGE', re.I)
# Whole-term matching correction 2026-11-23 (recheck scenario A): ^PRODUCT\b, POOR QUALITY (?!SLEEP), LACK OF(?! SPONTANEOUS SPEECH), (?<!SEVERE CUTANEOUS )ADVERSE REACTION$ so that
# Productive cough, Poor quality sleep (tesamorelin) and Lack of spontaneous speech, Severe cutaneous adverse reaction (somatropin comparator) are classed clinical.
# Audit correction 2026-11-23 (per protocol s7, dosing/medication-error terms): INCORRECT DOSAGE / DRUG ADMINISTRATION ERROR / DRUG TITRATION ERROR added to NONCLIN
# Post hoc (documented deviation): healthcare-utilisation/procedure terms reclassified as non-clinical after first output
POSTHOC_NONCLIN = {'EMERGENCY CARE','HOSPITALISATION','SURGERY','LIPOSUCTION'}
def nonclinical(pt, posthoc=True): return bool(NONCLIN.search(pt)) or (posthoc and pt.upper() in POSTHOC_NONCLIN)
# product/device/dose subset of non-clinical PTs (Objective e)
PRODUCT_DOSE = re.compile(r'PRODUCT|DEVICE|DOSE|MEDICATION ERROR|WRONG|INCORRECT|RECONSTITUTION|NEEDLE|SYRINGE|\bVIAL\b|DRUG STORAGE|EXPIRED|'
                          r'COUNTERFEIT|DISPENSING|PRESCRIBING|DRUG DELIVERY SYSTEM|WRONG TECHNIQUE|ACCIDENTAL|UNDERDOSE|OVERDOSE|INAPPROPRIATE SCHEDULE', re.I)
def product_dose(pt): return bool(PRODUCT_DOSE.search(pt))
# PTs that may reflect the treated condition / indication (flag only)
INDICATION_RELATED = {'LIPODYSTROPHY ACQUIRED','HIV LIPODYSTROPHY','LIPOHYPERTROPHY','ABDOMINAL DISTENSION','ABDOMINAL FAT INCREASED','FAT TISSUE INCREASED',
    'WEIGHT INCREASED','WEIGHT DECREASED','WAIST CIRCUMFERENCE INCREASED','OBESITY','HIV INFECTION','LIPOATROPHY','FAT REDISTRIBUTION','VISCERAL FAT INCREASED',
    'ABDOMINAL WALL DISORDER','CENTRAL OBESITY','BODY FAT DISORDER'}
# Post hoc additions (documented deviation): lipodystrophy/body-fat terms seen in first output
INDICATION_RELATED_POSTHOC = {'ABDOMINAL FAT APRON','FAT TISSUE DECREASED','SKIN LAXITY','FAT HYPERTROPHY','LIPOHYPERTROPHY'}
def indication_related(pt, posthoc=True): return pt.upper() in INDICATION_RELATED or (posthoc and pt.upper() in INDICATION_RELATED_POSTHOC)

# Labeled terms: Egrifta WR label (BLA 022505, effective 2026-07-29) sections 4-8; exact PTs plus documented synonyms
LABELED_EXACT = _exact(
  # 5.1 neoplasms handled by group; 5.2
  'Insulin-like growth factor increased','Blood growth hormone increased',
  # 5.3 fluid retention / 6.1 table
  'Oedema','Oedema peripheral','Peripheral swelling','Fluid retention','Generalised oedema','Swelling','Arthralgia','Pain in extremity','Myalgia',
  'Carpal tunnel syndrome','Joint swelling','Musculoskeletal stiffness','Musculoskeletal pain','Paraesthesia','Hypoaesthesia','Pain','Muscle strain',
  'Blood creatine phosphokinase increased','Night sweats','Palpitations','Dyspepsia','Vomiting',
  # 5.4 glucose
  'Hyperglycaemia','Blood glucose increased','Diabetes mellitus','Type 2 diabetes mellitus','Glucose tolerance impaired','Glycosylated haemoglobin increased',
  'Impaired fasting glucose','Diabetes mellitus inadequate control','Insulin resistance','Diabetic retinopathy','Retinopathy',
  # 5.5 hypersensitivity
  'Hypersensitivity','Drug hypersensitivity','Pruritus','Pruritus generalised','Erythema','Flushing','Urticaria','Rash','Rash pruritic','Rash erythematous',
  'Rash generalised','Rash macular','Rash papular',
  # 4/8 pregnancy; 7.2 glucocorticoid/adrenal
  'Exposure during pregnancy','Maternal exposure during pregnancy','Pregnancy','Adrenal insufficiency','Drug interaction')
LABELED_STARTS = ['INJECTION SITE']
def labeled(pt):
    p = pt.upper()
    if p in LABELED_EXACT or any(p.startswith(s) for s in LABELED_STARTS): return True
    return in_group(pt, 'Malignancy/neoplasm')  # label 5.1: increased risk of neoplasms
