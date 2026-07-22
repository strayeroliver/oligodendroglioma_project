# I got a lot of results with PAX5, so now I test if the methylation and rna appearances are correlated and whether methylation PAX5 holds up with the other variables independently
# Like age, has secondary oncogene, grade, etc. 

# this is to check if methylation data is actually valid to put inside the model otherwise its just riding the concordance of the RNA data

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import pandas as pd
import numpy as np
from lifelines import CoxPHFitter
from scipy import stats

oligo = pd.read_csv('data/TCGA/oligo_confirmed_master.csv')
oligo['grade_encoded'] = (oligo['neoplasm_histologic_grade'] == 'G3').astype(int)
oligo['secondary_og_burden'] = oligo[['NOTCH1', 'PIK3CA']].sum(axis=1)
oligo['has_secondary_og'] = (oligo['secondary_og_burden'] > 0).astype(int) # Just the same as before
rnaseq = pd.read_csv('data/TCGA/rnaseq_survival_ready.csv')
import json
with open('data/TCGA/key_gene_probes.json', 'r') as f:
    gene_probes = json.load(f)

methylation = pd.read_csv('data/TCGA/LGG_methylation450k.gz',
                           sep='\t', compression='gzip', index_col=0)

confirmed_ids = oligo['sampleID'].tolist()
ids_in_meth = [id for id in confirmed_ids if id in methylation.columns]
meth_oligo = methylation[ids_in_meth]

pax5_probes = [p for p in gene_probes['PAX5'] if p in meth_oligo.index]
print(f"PAX5 probes (in promoter): {len(pax5_probes)}")

pax5_meth = meth_oligo.loc[pax5_probes].mean(axis=0) * 100  
pax5_meth_df = pd.DataFrame({
    'sampleID': pax5_meth.index,
    'PAX5_methylation': pax5_meth.values
}) # Same process as before

merged = pax5_meth_df.merge(
    rnaseq[['sampleID', 'PAX5']].rename(columns={'PAX5': 'PAX5_expression'}),
    on='sampleID', how='inner'
).dropna()


r, p = stats.pearsonr(merged['PAX5_methylation'], merged['PAX5_expression'])
print(f"\nPAX5 methylation vs expression correlation:")
print(f"r = {r:.3f}, p = {p:.4f}") # These three were gemini snippits, I was confused on how to use scipy and stats command. 
# r is positive so they are correlated which does make sense. shared upstream regulatory mechanism.

# Scatter plot in matplotlib because I wanted to try it out (i only made bar graphs in my past with it before)
fig, ax = plt.subplots(figsize=(8, 6))
ax.scatter(merged['PAX5_methylation'], merged['PAX5_expression'], 
           alpha=0.6, color='#762a83')
ax.set_xlabel('PAX5 Promoter Methylation (%)', fontsize=12)
ax.set_ylabel('PAX5 Expression', fontsize=12)
ax.set_title(f'PAX5 Methylation vs Expression\nr={r:.3f}, p={p:.4f}', fontsize=12)
# My R value is 0.57 (R^2 = 0.325) but it makes sense as methylation data is definitely not the only thing that causes rna expression/limits it.
# Add trend line
z = np.polyfit(merged['PAX5_methylation'], merged['PAX5_expression'], 1) # https://numpy.org/doc/stable/reference/generated/numpy.polyfit.html
trend = np.poly1d(z)
x_range = np.linspace(merged['PAX5_methylation'].min(), 
                       merged['PAX5_methylation'].max(), 100)
ax.plot(x_range, trend(x_range), color='red', linestyle='--', linewidth=2)

plt.tight_layout()
plt.savefig('figures/08_PAX5_methylation_expression.png', dpi=300, bbox_inches='tight')



# Now we do full cox on methylation + clinical variables + oncogenes

full_data = oligo.merge(pax5_meth_df, on='sampleID', how='inner') # merge again for full dataset, we did all the hard work before
full_data = full_data[full_data['PFI.time'] > 0].dropna(
    subset=['PFI', 'PFI.time', 'PAX5_methylation']).copy()

full_data['PFI.time'] = full_data['PFI.time'] / 365 # years
multivariable_vars = ['PFI.time', 'PFI', 'grade_encoded', 
                       'age_at_initial_pathologic_diagnosis',
                       'has_secondary_og', 'PAX5_methylation'] 

cox_data = full_data[multivariable_vars].dropna()

cph = CoxPHFitter()
cph.fit(cox_data, duration_col='PFI.time', event_col='PFI')
print("\nFull test COX + Methylation data now")
cph.print_summary(decimals=3)

full_data.to_csv('data/TCGA/oligo_pax5_combined.csv', index=False)

# positively correlated so possible removal for model