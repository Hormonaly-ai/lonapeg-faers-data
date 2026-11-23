"""LONAPEG-FAERS figures from outputs/ (no new numbers)."""
import json, os, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); O = f'{HERE}/../outputs'; T = f'{O}/tables'; F = f'{O}/figures'; os.makedirs(F, exist_ok=True)
R = json.load(open(f'{O}/results.json')); V = json.load(open(f'{O}/review_results.json'))
plt.rcParams.update({'font.size': 9, 'font.family': 'DejaVu Sans'})
US = {'Hyperglycaemia/type 2 diabetes': 'Hyperglycemia/type 2 diabetes', 'CNS/pituitary tumour': 'CNS/pituitary tumor', 'NC: any non-clinical': 'Non-clinical: any', 'NC: Dose omission/missed or delayed dose': 'Non-clinical: dose omission/missed or delayed dose', 'NC: Device issues': 'Non-clinical: device issues', 'NC: Storage/temperature/expiry/preparation': 'Non-clinical: storage/temperature/expiry/preparation', 'NC: Incorrect dose/technique/schedule/product': 'Non-clinical: incorrect dose/technique/schedule/product', 'NC: Product quality/supply/access': 'Non-clinical: product quality/supply/access'}
us = lambda s: US.get(s, s)
# Fig 1: quarterly reports
q = pd.read_csv(f'{T}/table_quarterly_counts.csv')
fig, ax = plt.subplots(figsize=(8, 4.2))
x = np.arange(len(q)); b1 = q.L_PS_ped; b2 = q.L_PS_adult; b3 = q.L_PS_age_unknown
ax.bar(x, b1, color='#2b6cb0', label='Age <18 years'); ax.bar(x, b2, bottom=b1, color='#c05621', label='Age ≥18 years'); ax.bar(x, b3, bottom=b1 + b2, color='#a0aec0', label='Age not reported')
ax.plot(x, q.L_PS_only_nonclinical_PTs, 'k--o', ms=3, lw=1, label='Reports with only product/device/dosing PTs')
ax.set_xticks(x); ax.set_xticklabels(q.quarter, rotation=60, ha='right', fontsize=7.5); ax.set_ylabel('Primary-suspect reports (n)')
ax.axvline(list(q.quarter).index('2025Q3') - 0.5, color='#c05621', lw=0.8, ls=':'); ax.text(list(q.quarter).index('2025Q3') - 0.4, ax.get_ylim()[1] * 0.92, 'US adult\napproval\n(Jul 2025)', fontsize=7, color='#c05621')
ax.legend(frameon=False, fontsize=7.5, loc='upper left'); ax.spines[['top', 'right']].set_visible(False)
plt.tight_layout(); plt.savefig(f'{F}/Fig1.png', dpi=300); plt.close()
# Fig 2: forest, prespecified clinical groups + NC groups vs all drugs
g = pd.read_csv(f'{T}/table3_prespecified_groups_L_PS.csv')
g = g[g.a >= 1].iloc[::-1].reset_index(drop=True)
fig, ax = plt.subplots(figsize=(7.2, 0.32 * len(g) + 1))
for i, r in g.iterrows():
    c = '#c53030' if r.signal_primary else '#4a5568'
    ax.plot([r.ROR_lo, r.ROR_hi], [i, i], color=c, lw=1.2); ax.plot(r.ROR, i, 's', color=c, ms=4)
    ax.text(600, i, f"{int(r.a)}  {r.ROR:.2f} ({r.ROR_lo:.2f}–{r.ROR_hi:.2f})", va='center', fontsize=7.5)
ax.axvline(1, color='k', lw=0.7); ax.set_xscale('log'); ax.set_xlim(0.01, 500)
ax.set_yticks(range(len(g))); ax.set_yticklabels([us(x) for x in g.group], fontsize=7.5)
ax.set_xlabel('Reporting odds ratio vs all other reports (95% CI, log scale)'); ax.text(600, len(g) - 0.2, 'n   ROR (95% CI)', fontsize=7.5, fontweight='bold')
ax.spines[['top', 'right']].set_visible(False); plt.tight_layout(); plt.savefig(f'{F}/Fig2.png', dpi=300, bbox_inches='tight'); plt.close()
# Fig 3: comparators (primary and S12) for clinical groups with any count
c = pd.read_csv(f'{T}/table4_active_comparators.csv'); s = pd.read_csv(f'{T}/tableS_S12_comparators_clinical_reports_only.csv')
items = ['Injection-site reactions', 'Lipoatrophy/lipohypertrophy', 'Hypersensitivity', 'Neoplasm (any)', 'Hyperglycaemia/type 2 diabetes', 'IGF-1 increased', 'Scoliosis', 'Intracranial hypertension', 'Fluid retention/musculoskeletal', 'Headache']
fig, axes = plt.subplots(1, 2, figsize=(9, 4.8), sharey=True)
for ax, (an_p, an_s, title) in zip(axes, [('L_PS vs DS_PS', 'S12 (post hoc): L_PS vs DS_PS, reports with >=1 clinical PT', 'vs daily somatropin'), ('L_PS vs SG_PS', 'S12 (post hoc): L_PS vs SG_PS, reports with >=1 clinical PT', 'vs somatrogon')]):
    for j, it in enumerate(items[::-1]):
        for off, d, col in [(0.15, c[c.analysis == an_p], '#2b6cb0'), (-0.15, s[s.analysis == an_s], '#dd6b20')]:
            r = d[d['item'] == it].iloc[0]
            ax.plot([max(r.ROR_lo, 0.005), min(r.ROR_hi, 900)], [j + off] * 2, color=col, lw=1.1); ax.plot(r.ROR, j + off, 's' if r.c_comp > 0 and r.a_L > 0 else 'D', color=col, ms=3.5, mfc=col if r.c_comp > 0 and r.a_L > 0 else 'white')
    ax.axvline(1, color='k', lw=0.7); ax.set_xscale('log'); ax.set_xlim(0.01, 1000); ax.set_title(f'Lonapegsomatropin {title}', fontsize=9)
    ax.set_xlabel('ROR (95% CI, log scale)'); ax.spines[['top', 'right']].set_visible(False)
axes[0].set_yticks(range(len(items))); axes[0].set_yticklabels([us(x) for x in items[::-1]], fontsize=7.5)
h1, = axes[1].plot([], [], color='#2b6cb0', marker='s', label='All primary-suspect reports (prespecified)'); h2, = axes[1].plot([], [], color='#dd6b20', marker='s', label='Reports with ≥1 clinical PT (post hoc S12)')
h3, = axes[1].plot([], [], color='grey', marker='D', mfc='white', ls='none', label='Zero comparator reports (Haldane correction)')
fig.legend(handles=[h1, h2, h3], frameon=False, fontsize=7.5, loc='lower center', ncol=3); plt.tight_layout(rect=(0, 0.07, 1, 1)); plt.savefig(f'{F}/Fig3.png', dpi=300); plt.close()
# Fig 4: PRISMA flow
P = V['prisma']
fig, ax = plt.subplots(figsize=(8, 7)); ax.axis('off')
def box(x, y, w, h, t):
    ax.add_patch(plt.Rectangle((x, y), w, h, fill=False, lw=0.9)); ax.text(x + w / 2, y + h / 2, t, ha='center', va='center', fontsize=7.6, wrap=True)
def arr(x1, y1, x2, y2): ax.annotate('', xy=(x2, y2), xytext=(x1, y1), arrowprops=dict(arrowstyle='->', lw=0.8))
box(0.05, 0.84, 0.42, 0.12, f"Records identified (search {P['search_date']}):\nPubMed n = {P['records_pubmed']}\nEurope PMC (incl. abstracts, preprints) n = {P['records_epmc']}")
box(0.55, 0.84, 0.40, 0.12, f"Duplicates removed n = {P['duplicates_removed']}\n\nClinicalTrials.gov: {P['ctgov_records']} registrations\nchecked; no additional studies")
box(0.05, 0.66, 0.42, 0.10, f"Records screened (title/abstract)\nn = {P['records_screened']}"); box(0.55, 0.66, 0.40, 0.10, f"Records excluded\nn = {P['records_excluded_title_abstract']}")
box(0.05, 0.50, 0.42, 0.10, f"Reports sought for retrieval\nn = {P['reports_sought']}"); box(0.55, 0.50, 0.40, 0.10, f"Reports not retrieved\nn = {P['reports_not_retrieved']} (no institutional license)")
ex = '\n'.join(f"{v}: {k.split(' (')[0]}" for k, v in P['excluded_reasons'].items())
box(0.05, 0.32, 0.42, 0.12, f"Reports assessed for eligibility (full text)\nn = {P['reports_assessed']}"); box(0.55, 0.28, 0.40, 0.18, f"Reports excluded n = {P['reports_excluded']}\n{ex}")
box(0.05, 0.08, 0.42, 0.16, f"Reports included n = {P['reports_included']}\nStream 1 (cases/series/cohort) n = {P['stream1_reports']}\nStream 2 (trials) n = {P['stream2_reports']} ({P['stream2_unique_studies']} studies)\n+ regulatory documents n = {P['regulatory_documents']}")
for y1, y2 in [(0.84, 0.76), (0.66, 0.60), (0.50, 0.44), (0.32, 0.24)]: arr(0.26, y1, 0.26, y2)
for y in [0.90, 0.71, 0.55, 0.38]: arr(0.47, y, 0.55, y)
plt.savefig(f'{F}/Fig4.png', dpi=300, bbox_inches='tight'); plt.close()
print('figures ok')
