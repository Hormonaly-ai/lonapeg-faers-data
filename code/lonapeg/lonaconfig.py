"""Prespecified definitions for LONAPEG-FAERS (protocol v1.0, OSF 10.17605/OSF.IO/UB6W2). Adapted from tesaconfig.py (TESA-FAERS)."""
import re
WINDOW = [f'{y}Q{q}' for y in range(2021, 2027) for q in range(1, 5) if (y, q) >= (2021, 3) and (y, q) <= (2026, 3)]
LONA = re.compile(r'LONAPEGSOMATROPIN|SKYTROFA|TRANSCON\s*(H?GH|GROWTH)|ACP-?011', re.I)
SOMATROGON = re.compile(r'SOMATROGON|NGENLA', re.I)
SOMAPACITAN = re.compile(r'SOMAPACITAN|SOGROYA', re.I)
SOM_ANY = re.compile(r'SOMATROPIN|GENOTROPIN|HUMATROPE|NORDITROPIN|NUTROPIN|OMNITROPE|SAIZEN|ZOMACTON|TEV-TROPIN|ACCRETROPIN|VALTROPIN|SEROSTIM|ZORBTIVE', re.I)
SERO_ZORB = re.compile(r'SEROSTIM|ZORBTIVE', re.I)
GH_ALL = re.compile(LONA.pattern + '|' + SOMATROGON.pattern + '|' + SOMAPACITAN.pattern + '|' + SOM_ANY.pattern, re.I)
def drug_class(k):
    k = k or ''
    if LONA.search(k): return 'LONA'
    if SOMATROGON.search(k): return 'SOMATROGON'
    if SOMAPACITAN.search(k): return 'SOMAPACITAN'
    if SERO_ZORB.search(k): return 'SERO_ZORB'
    if SOM_ANY.search(k): return 'DAILY_SOM'
    return None

def _exact(*xs): return {x.upper() for x in xs}
NEO_ANY = ['NEOPLASM','TUMOUR','CARCINOMA','CANCER','LYMPHOMA','LEUKAEMIA','SARCOMA','MELANOMA','METASTA','MYELOMA','MALIGNAN','ADENOMA','MENINGIOMA',
           'CRANIOPHARYNGIOMA','GLIOMA','ASTROCYTOMA','MEDULLOBLASTOMA','GERMINOMA','EPENDYMOMA','BLASTOMA']
NEO_EXCLUDE = re.compile(r'TUMOUR LYSIS|TUMOUR MARKER|NEOPLASM PROPHYLAXIS|BLASTOMYCOSIS|TUMOUR PAIN', re.I)
MALIG = ['MALIGNAN','CARCINOMA','CANCER','LYMPHOMA','LEUKAEMIA','SARCOMA','MELANOMA','METASTA','MYELOMA']
CNS = ['BRAIN','PITUITARY','CRANIOPHARYNGIOMA','GLIOMA','MENINGIOMA','ASTROCYTOMA','MEDULLOBLASTOMA','GERMINOMA','EPENDYMOMA','CENTRAL NERVOUS SYSTEM']
GROUPS = {
 'Injection-site reactions': dict(startswith=['INJECTION SITE']),
 'Lipoatrophy/lipohypertrophy': dict(exact=_exact('Lipoatrophy','Injection site atrophy','Lipohypertrophy','Lipodystrophy acquired','Fat tissue decreased',
    'Skin atrophy','Injection site lipoatrophy','Injection site lipohypertrophy')),
 'Hypersensitivity': dict(exact=_exact('Hypersensitivity','Drug hypersensitivity','Anaphylactic reaction','Anaphylactic shock','Anaphylactoid reaction',
    'Urticaria','Rash','Rash pruritic','Rash erythematous','Rash generalised','Rash macular','Rash papular','Pruritus','Pruritus generalised',
    'Erythema','Angioedema','Swelling face','Lip swelling','Flushing')),
 'Neoplasm (any)': dict(neo_any=True),
 'Malignancy': dict(contains=MALIG, neoplasm_not_benign=True),
 'CNS/pituitary tumour': dict(cns=True),
 'Hyperglycaemia/type 2 diabetes': dict(exact=_exact('Hyperglycaemia','Blood glucose increased','Diabetes mellitus','Type 2 diabetes mellitus',
    'Diabetes mellitus inadequate control','Glycosylated haemoglobin increased','Glucose tolerance impaired','Impaired fasting glucose','Insulin resistance',
    'Diabetic ketoacidosis','Blood glucose abnormal')),
 'IGF-1 increased': dict(exact=_exact('Insulin-like growth factor increased','Blood growth hormone increased')),
 'SCFE/osteonecrosis': dict(exact=_exact('Epiphysiolysis','Slipped capital femoral epiphysis','Osteonecrosis','Osteonecrosis of femoral head','Perthes disease','Avascular necrosis')),
 'Scoliosis': dict(exact=_exact('Scoliosis','Kyphoscoliosis','Kyphosis','Spinal deformity')),
 'Intracranial hypertension': dict(exact=_exact('Idiopathic intracranial hypertension','Intracranial pressure increased','Benign intracranial hypertension','Papilloedema')),
 'Fluid retention/musculoskeletal': dict(exact=_exact('Oedema peripheral','Oedema','Peripheral swelling','Fluid retention','Swelling','Generalised oedema',
    'Arthralgia','Pain in extremity','Myalgia','Joint swelling','Musculoskeletal stiffness','Musculoskeletal pain','Paraesthesia','Hypoaesthesia','Carpal tunnel syndrome')),
 'Hypothyroidism': dict(exact=_exact('Hypothyroidism','Central hypothyroidism','Secondary hypothyroidism','Thyroxine free decreased','Thyroxine decreased',
    'Blood thyroid stimulating hormone decreased','Tri-iodothyronine free decreased')),
 'Adrenal insufficiency': dict(exact=_exact('Adrenal insufficiency','Secondary adrenocortical insufficiency','Adrenocortical insufficiency acute','Blood cortisol decreased','Glucocorticoid deficiency')),
 'Pancreatitis': dict(exact=_exact('Pancreatitis','Pancreatitis acute','Pancreatitis chronic')),
 'Headache': dict(exact=_exact('Headache')),
 'Pregnancy exposure': dict(exact=_exact('Exposure during pregnancy','Maternal exposure during pregnancy','Pregnancy','Abortion spontaneous')),
}
def in_group(pt, g):
    p = pt.upper(); d = GROUPS[g]
    if p in d.get('exact', set()): return True
    if any(p.startswith(s) for s in d.get('startswith', [])): return True
    if d.get('neo_any'): return any(s in p for s in NEO_ANY) and not NEO_EXCLUDE.search(p)
    if d.get('cns'): return any(s in p for s in NEO_ANY) and not NEO_EXCLUDE.search(p) and any(s in p for s in CNS)
    if any(s in p for s in d.get('contains', [])): return True
    if d.get('neoplasm_not_benign') and 'NEOPLASM' in p and 'BENIGN' not in p and not NEO_EXCLUDE.search(p): return True
    return False

# Non-clinical PTs (protocol s7): reused from TESA-FAERS tesaconfig.NONCLIN (incl. its 2026-11-23 audit correction), plus device/cartridge/storage terms
NONCLIN = re.compile(r'OFF LABEL|PRODUCT |^PRODUCT|PRODUCT$|DRUG INEFFECTIVE|WRONG TECHNIQUE|NO ADVERSE EVENT|INTENTIONAL PRODUCT|DEVICE|DOSE OMISSION|'
                     r'INAPPROPRIATE SCHEDULE|INCORRECT DOSE|EXTRA DOSE|WRONG (DRUG|PATIENT|PRODUCT|DOSE|ROUTE)|MEDICATION ERROR|EXPIRED PRODUCT|'
                     r'THERAPEUTIC RESPONSE (UNEXPECTED|DECREASED)|DRUG EFFECT (LESS|DECREASED|DELAYED|INCOMPLETE)|TREATMENT NONCOMPLIANCE|NEEDLE ISSUE|'
                     r'UNDERDOSE|ACCIDENTAL EXPOSURE|CIRCUMSTANCE OR INFORMATION|ADVERSE EVENT$|^ILLNESS$|INSURANCE|PRESCRIBED|DRUG DOSE OMISSION|'
                     r'INTERCEPTED|SUSPECTED COUNTERFEIT|DRUG DIVERSION|NO REACTION|UNEVALUABLE EVENT|CONDITION AGGRAVATED|DRUG ABUSE|'
                     r'SYRINGE ISSUE|INJURY ASSOCIATED WITH DEVICE|LIQUID PRODUCT|DRUG ADMINISTERED|LACK OF|DRUG DISPENSING|UNINTENTIONAL|'
                     r'RECONSTITUTION|POOR QUALITY|DRUG DOSE TITRATION|DOSE CALCULATION|DRUG PRESCRIBING|INCORRECT ROUTE|'
                     r'INCORRECT PRODUCT|MISSED DOSE|INCORRECT (DRUG )?DOSAGE|ADMINISTRATION ERROR|TITRATION ERROR|DRUG DELIVERY SYSTEM|ADVERSE DRUG REACTION$|ADVERSE REACTION$|THERAPY INTERRUPTED|THERAPY CESSATION|'
                     r'THERAPY NON-RESPONDER|TREATMENT FAILURE|DRUG INTOLERANCE|INTENTIONAL DOSE|OVERDOSE|\bVIAL\b|DRUG STORAGE|'
                     r'CARTRIDGE|INJECTOR|INCORRECT DRUG ADMINISTRATION|DRUG PREPARATION|HANDLING ERROR|STORAGE|TEMPERATURE EXCURSION|'
                     r'DRUG ADMINISTRATION ERROR|INAPPROPRIATE|DELAYED DOSE|DRUG DOSING|DOSE NOT|INCORRECT FREQUENCY|DRUG SUPPLY|INTERRUPTION|'
                     r'TREATMENT DELAYED|THERAPY CHANGE|DRUG LEVEL|ACCIDENTAL UNDERDOSE|ACCIDENTAL OVERDOSE', re.I)
# PTs not matching NONCLIN but non-clinical by nature (procedural/administrative), reviewed at first output (any additions = documented deviation)
POSTHOC_NONCLIN = set()
# Correction (2026-11-23, before final analysis; reported as a deviation): clinical PTs that the inherited text pattern mis-captured
CLIN_EXCEPT = {'PRODUCTIVE COUGH','POOR QUALITY SLEEP','INAPPROPRIATE ANTIDIURETIC HORMONE SECRETION','INAPPROPRIATE AFFECT','INAPPROPRIATE SINUS TACHYCARDIA'}
def nonclinical(pt, posthoc=True):
    if pt.upper() in CLIN_EXCEPT: return False
    return bool(NONCLIN.search(pt)) or (posthoc and pt.upper() in POSTHOC_NONCLIN)
NC_SUB = [
 ('Dose omission/missed or delayed dose', re.compile(r'DOSE OMISSION|MISSED DOSE|DELAYED DOSE|DOSE NOT|THERAPY INTERRUPTED|INTERRUPTION|TREATMENT NONCOMPLIANCE|TREATMENT DELAYED|THERAPY CESSATION', re.I)),
 ('Device issues', re.compile(r'DEVICE|INJECTOR|CARTRIDGE|NEEDLE|SYRINGE|DRUG DELIVERY SYSTEM', re.I)),
 ('Storage/temperature/expiry/preparation', re.compile(r'STORAGE|TEMPERATURE|EXPIRED|RECONSTITUTION|PREPARATION|HANDLING|LIQUID PRODUCT|PRODUCT PHYSICAL|PRODUCT COLOUR|PRODUCT CONTAMINATION|PRODUCT LEAKAGE', re.I)),
 ('Incorrect dose/technique/schedule/product', re.compile(r'INCORRECT|WRONG|EXTRA DOSE|UNDERDOSE|OVERDOSE|INAPPROPRIATE SCHEDULE|INCORRECT FREQUENCY|ADMINISTRATION ERROR|MEDICATION ERROR|ACCIDENTAL|DOSE CALCULATION|TITRATION|DRUG DOSING|INTENTIONAL', re.I)),
 ('Product quality/supply/access', re.compile(r'QUALITY|COUNTERFEIT|DRUG SUPPLY|SUPPLY ISSUE|AVAILABILITY|INSURANCE|PRESCRIBED|PRODUCT DISTRIBUTION|DRUG DISPENSING', re.I)),
]
def nc_subcats(pt): return [n for n, r in NC_SUB if r.search(pt)]
INDICATION_RELATED = {'GROWTH RETARDATION','SHORT STATURE','BODY HEIGHT BELOW NORMAL','GROWTH HORMONE DEFICIENCY','INSULIN-LIKE GROWTH FACTOR DECREASED',
    'BLOOD GROWTH HORMONE DECREASED','GROWTH FAILURE','HYPOPITUITARISM','DELAYED PUBERTY','BODY HEIGHT ABNORMAL','GROWTH ACCELERATED','BODY HEIGHT INCREASED'}
INDICATION_RELATED_POSTHOC = set()
def indication_related(pt, posthoc=True): return pt.upper() in INDICATION_RELATED or (posthoc and pt.upper() in INDICATION_RELATED_POSTHOC)

# Labeled terms: Skytrofa US label (BLA 761177, openFDA version 14, effective 2026-10-21), sections 4-8 (protocol s11)
LABELED_EXACT = _exact(
  'Hypersensitivity','Drug hypersensitivity','Anaphylactic reaction','Anaphylactic shock','Anaphylactoid reaction','Angioedema','Urticaria','Rash',
  'Rash pruritic','Rash erythematous','Rash generalised','Rash macular','Rash papular','Pruritus','Pruritus generalised','Erythema','Flushing','Swelling face','Lip swelling',
  'Hyperglycaemia','Blood glucose increased','Diabetes mellitus','Type 2 diabetes mellitus','Diabetes mellitus inadequate control','Glycosylated haemoglobin increased',
  'Glucose tolerance impaired','Impaired fasting glucose','Insulin resistance','Diabetic ketoacidosis','Blood glucose abnormal','Diabetic retinopathy',
  'Idiopathic intracranial hypertension','Intracranial pressure increased','Benign intracranial hypertension','Papilloedema','Visual impairment','Vision blurred','Headache','Nausea','Vomiting',
  'Oedema','Oedema peripheral','Peripheral swelling','Fluid retention','Generalised oedema','Swelling','Arthralgia','Myalgia','Carpal tunnel syndrome','Paraesthesia','Arthritis','Arthritis reactive',
  'Adrenal insufficiency','Secondary adrenocortical insufficiency','Blood cortisol decreased','Hypothyroidism','Central hypothyroidism','Secondary hypothyroidism',
  'Thyroxine free decreased','Thyroxine decreased','Blood thyroid stimulating hormone decreased','Tri-iodothyronine free decreased',
  'Epiphysiolysis','Slipped capital femoral epiphysis','Osteonecrosis','Osteonecrosis of femoral head','Gait disturbance','Limp','Scoliosis','Kyphoscoliosis',
  'Pancreatitis','Pancreatitis acute','Abdominal pain','Abdominal pain upper','Lipoatrophy','Sudden death',
  'Viral infection','Pyrexia','Cough','Haemorrhage','Epistaxis','Contusion','Petechiae','Eye haemorrhage','Diarrhoea',
  'Blood alkaline phosphatase increased','Blood phosphorus increased','Blood parathyroid hormone increased',
  'Insulin-like growth factor increased','Naevus','Melanocytic naevus','Meningioma')
LABELED_STARTS = ['INJECTION SITE']
def labeled(pt):
    p = pt.upper()
    if p in LABELED_EXACT or any(p.startswith(s) for s in LABELED_STARTS): return True
    return in_group(pt, 'Neoplasm (any)')  # label 5.3 neoplasms
