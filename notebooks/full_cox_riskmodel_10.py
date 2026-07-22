# Combined multi-layer risk model
# Integrating clinical, genomic, transcriptomic, and epigenetic predictors into a single  score

# This whole file is really just copy and pasting all the final fragments from every file while doing cox repeatedly.

# This is the foundation for the clinical prediction tool (nomogram) that I talked w mentor

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import pandas as pd
import numpy as np
import json
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
from scipy import stats



# Clinical + genomic 
oligo = pd.read_csv('data/TCGA/oligo_confirmed_master.csv')
oligo['grade_encoded'] = (oligo['neoplasm_histologic_grade'] == 'G3').astype(int)
oligo['secondary_og_burden'] = oligo[['NOTCH1', 'PIK3CA']].sum(axis=1)
oligo['has_secondary_og'] = (oligo['secondary_og_burden'] > 0).astype(int)
rnaseq = pd.read_csv('data/TCGA/rnaseq_survival_ready.csv')
lasso_coefs = {
    'VEPH1': 0.18971, 'ADAM6': 0.11343, 'TRH': 0.20784,
    'HLA-DQA2': 0.11264, 'SIX1': 0.04329, 'CENPV': -0.26925,
    'ABCC3': 0.25340, 'DLX6': 0.22175, 'C18orf34': -0.27209,
    'PAX5': 0.36958, 'SEL1L3': 0.11913
}
rnaseq['rnaseq_risk_score'] = 0.0
for gene, beta in lasso_coefs.items():
    if gene in rnaseq.columns:
        rnaseq['rnaseq_risk_score'] += rnaseq[gene] * beta

# PAX5 methylation (promoter from notebook 09)
with open('data/TCGA/key_gene_probes.json', 'r') as f:
    gene_probes = json.load(f)
methylation = pd.read_csv('data/TCGA/LGG_methylation450k.gz',
                           sep='\t', compression='gzip', index_col=0)
confirmed_ids = oligo['sampleID'].tolist()
ids_in_meth = [id for id in confirmed_ids if id in methylation.columns]
meth_oligo = methylation[ids_in_meth]
pax5_probes = [p for p in gene_probes['PAX5'] if p in meth_oligo.index]
pax5_meth = meth_oligo.loc[pax5_probes].mean(axis=0) * 100
pax5_df = pd.DataFrame({'sampleID': pax5_meth.index, 
                         'PAX5_methylation': pax5_meth.values})

# Merge into a dataset
combined = oligo.merge(
    rnaseq[['sampleID', 'rnaseq_risk_score']], on='sampleID', how='inner'
).merge(
    pax5_df, on='sampleID', how='inner'
)

combined = combined[combined['PFI.time'] > 0].dropna(
    subset=['PFI', 'PFI.time', 'grade_encoded', 
            'age_at_initial_pathologic_diagnosis',
            'has_secondary_og', 'rnaseq_risk_score', 'PAX5_methylation']
).copy()
combined['PFI.time'] = combined['PFI.time'] / 365


# cox it up
model_vars = ['PFI.time', 'PFI', 'grade_encoded', 
              'age_at_initial_pathologic_diagnosis',
              'has_secondary_og', 'rnaseq_risk_score', 
              'PAX5_methylation']

cox_data = combined[model_vars].dropna()
print(f"Complete cases for model: {len(cox_data)}")

cph_combined = CoxPHFitter(penalizer=0.1)
cph_combined.fit(cox_data, duration_col='PFI.time', event_col='PFI')

print("\nFull model w/ everything")
cph_combined.print_summary(decimals=3)

print("\nModel concordance (we add 1 additional variable every time, should see upwards trend)")

# Clinical only
cph_clin = CoxPHFitter()
cph_clin.fit(cox_data[['PFI.time', 'PFI', 'grade_encoded', 
                        'age_at_initial_pathologic_diagnosis']],
             duration_col='PFI.time', event_col='PFI')
print(f"Clinical only:                    {cph_clin.concordance_index_:.3f}") # Rediculous # of spaces to make it look nice and pretty in terminal

# Clinical + genomic
cph_genomic = CoxPHFitter()
cph_genomic.fit(cox_data[['PFI.time', 'PFI', 'grade_encoded',
                           'age_at_initial_pathologic_diagnosis', 
                           'has_secondary_og']],
                duration_col='PFI.time', event_col='PFI')
print(f"Clinical + genomic:               {cph_genomic.concordance_index_:.3f}")

# Clinical + genomic + transcriptomic
cph_transcript = CoxPHFitter(penalizer=0.1)
cph_transcript.fit(cox_data[['PFI.time', 'PFI', 'grade_encoded',
                              'age_at_initial_pathologic_diagnosis',
                              'has_secondary_og', 'rnaseq_risk_score']],
                    duration_col='PFI.time', event_col='PFI')
print(f"Clinical + genomic + transcript:  {cph_transcript.concordance_index_:.3f}")

# All four layers
print(f"All four layers combined:         {cph_combined.concordance_index_:.3f}")

combined['final_risk_score'] = cph_combined.predict_partial_hazard(
    cox_data[['grade_encoded', 'age_at_initial_pathologic_diagnosis',
              'has_secondary_og', 'rnaseq_risk_score', 'PAX5_methylation']]
) # Done before, we see an increase of concordance every single time it runs which is really good.

median_risk = combined['final_risk_score'].median()
combined['final_risk_group'] = (combined['final_risk_score'] > median_risk).astype(int)

high = combined[combined['final_risk_group'] == 1]
low = combined[combined['final_risk_group'] == 0]

lr = logrank_test(
    high['PFI.time'], low['PFI.time'],
    event_observed_A=high['PFI'], event_observed_B=low['PFI']
)

fig, ax = plt.subplots(figsize=(10, 6))
kmf_h = KaplanMeierFitter()
kmf_h.fit(high['PFI.time'], event_observed=high['PFI'],
          label=f'High Risk (n={len(high)})')
kmf_h.plot_survival_function(ax=ax, ci_show=True, color='#d7191c')

kmf_l = KaplanMeierFitter()
kmf_l.fit(low['PFI.time'], event_observed=low['PFI'],
          label=f'Low Risk (n={len(low)})')
kmf_l.plot_survival_function(ax=ax, ci_show=True, color='#2c7bb6')

ax.set_xlabel('Time (Years)', fontsize=12)
ax.set_ylabel('Progression-Free Survival', fontsize=12)
ax.set_title(f'Multi-Layer Risk Score\nClinical + Genomic + Transcriptomic + Epigenetic',
             fontsize=12)
plt.tight_layout()
plt.savefig('figures/10_final_multilayer_risk_score.png', dpi=300, bbox_inches='tight')

print(f"High risk: {kmf_h.median_survival_time_:.2f} years")
print(f"Low risk: {kmf_l.median_survival_time_:.2f} years")

combined.to_csv('data/TCGA/oligo_final_combined_model.csv', index=False) 


