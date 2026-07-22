# This file will be about labeling the genes as tumor suppresors vs secondary oncogenes and looking at the ratios between them 
# Should be pretty short, just using the OncoKB database 

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import pandas as pd
import numpy as np
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines.statistics import proportional_hazard_test


oligo = pd.read_csv('data/TCGA/oligo_confirmed_master.csv') # Load in confirmed dataset validated by cBioPortal from cox_regression_02

# Tumor suppressors are loss of function mutations
# They suppress growth, their loss essentially removes the brakes of the tumor growth (so a loss of them means you have increased risk)

# Secondary Oncogenes are gain of function mutations
# They are essentially the gas pedal, more oncogenes and less suppressors mean even higher risk

tumor_suppressors = ['CIC', 'FUBP1', 'ATRX'] # https://www.oncokb.org/ for all classification

secondary_oncogenes = ['NOTCH1', 'PIK3CA'] # I removed IDH 1/2 because IDH is in 100% of all oligodendroglioma

available_ts = [g for g in tumor_suppressors if g in oligo.columns]
available_og = [g for g in secondary_oncogenes if g in oligo.columns]

oligo['ts_burden'] = oligo[available_ts].sum(axis=1) # Summing up all the tumor suppressors/oncogenes to get a 'burden' value (how many there are)
oligo['og_burden'] = oligo[available_og].sum(axis=1) # higher burden value - more of that type
oligo['has_secondary_og'] = (oligo['og_burden'] > 0).astype(int) # binary has or doesnt have
oligo['has_ts_mutation'] = (oligo['ts_burden'] > 0).astype(int)

print(f"Mean TS: {oligo['ts_burden'].mean():.3f}") # We see 0.953 TS burden, 0.339 OG burden. Makes sense, there are traditionally always more TS > OG
print(f"Mean OG: {oligo['og_burden'].mean():.3f}") # 37 w OG, 93 w TS
print(f"Patients with OG mutation: {oligo['has_secondary_og'].sum()}")
print(f"Patients with TS mutation: {oligo['has_ts_mutation'].sum()}")

# Filtering PFI data like always
pfi_data = oligo[oligo['PFI.time'] > 0].dropna(
    subset=['PFI', 'PFI.time']
).copy()
pfi_data['PFI.time'] = pfi_data['PFI.time'] / 365  # years this time, I'm makign all my graphs be in years

# Cox:
cox_variables = {
    'Secondary OG mutation (any)': 'has_secondary_og',
    'TS mutation (any)': 'has_ts_mutation',
    'Secondary OG burden (count)': 'og_burden',
    'TS burden (count)': 'ts_burden'
}
cph = CoxPHFitter()
results = []

for name, col in cox_variables.items():
    data = pfi_data[['PFI.time', 'PFI', col]].dropna()
    if data[col].nunique() < 2: # skipping again if lacks varaince
        continue

    cph.fit(data, duration_col='PFI.time', event_col='PFI')
    s = cph.summary # Exact same format in CGGA_validation
    hr = s['exp(coef)'].values[0]
    ci_low = s['exp(coef) lower 95%'].values[0]
    ci_high = s['exp(coef) upper 95%'].values[0]
    p = s['p'].values[0]
    results.append({
        'Variable': name,
        'HR': round(hr, 3),
        '95% CI': f"{round(ci_low, 2)}-{round(ci_high, 2)}", 
        'p-value': round(p, 4)
    })

results_df = pd.DataFrame(results)
print(results_df.to_string(index=False))
results_df.to_csv('results/oncokb_cox_results.csv', index=False)



# Now make KM graphs to show the differences visually in survival between those with and without

og_pos = pfi_data[pfi_data['has_secondary_og'] == 1]
og_neg = pfi_data[pfi_data['has_secondary_og'] == 0]

lr = logrank_test( # Same log rank p test as before: https://numiqo.com/tutorial/log-rank-test
    og_pos['PFI.time'], og_neg['PFI.time'],
    event_observed_A=og_pos['PFI'],
    event_observed_B=og_neg['PFI']
)

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# We will be plotting two lines 
kmf_pos = KaplanMeierFitter()
kmf_pos.fit(og_pos['PFI.time'], event_observed=og_pos['PFI'],
            label=f'Secondary OG mutant (n={len(og_pos)})')
kmf_pos.plot_survival_function(ax=axes[0], ci_show=True, color='#d7191c')

kmf_neg = KaplanMeierFitter()
kmf_neg.fit(og_neg['PFI.time'], event_observed=og_neg['PFI'],
            label=f'No secondary OG (n={len(og_neg)})')
kmf_neg.plot_survival_function(ax=axes[0], ci_show=True, color='#2c7bb6')

axes[0].set_xlabel('Time (Years)', fontsize=12)
axes[0].set_ylabel('Progression-Free Survival', fontsize=12)
axes[0].set_title('Secondary Oncogene Mutation vs PFI\nMolecularly Confirmed Oligodendroglioma', fontsize=11)
axes[0].text(0.02, 0.05, f'Log-rank p={lr.p_value:.4f}',
             transform=axes[0].transAxes, fontsize=10)
axes[0].grid(True, linestyle=':', alpha=0.5)

print(f"\nSecondary OG median PFI: {kmf_pos.median_survival_time_:.2f} years")
print(f"No secondary OG median PFI: {kmf_neg.median_survival_time_:.2f} years")
print(f"Log-rank p: {lr.p_value:.4f}") # https://numiqo.com/tutorial/log-rank-test

# This graph will be for tumor suppressors
pfi_data['ts_group'] = pfi_data['ts_burden'].clip(upper=2)
colors_ts = ['#2c7bb6', '#fdae61', '#d7191c']
labels_ts = ['0 TS mutations', '1 TS mutation', '2+ TS mutations']

for group, color, label in zip([0, 1, 2], colors_ts, labels_ts):
    subset = pfi_data[pfi_data['ts_group'] == group]
    if len(subset) < 5:
        continue
    kmf = KaplanMeierFitter()
    kmf.fit(subset['PFI.time'], event_observed=subset['PFI'],
            label=f'{label} (n={len(subset)})')
    kmf.plot_survival_function(ax=axes[1], ci_show=False, color=color)

axes[1].set_xlabel('Time (Years)', fontsize=12)
axes[1].set_ylabel('Progression-Free Survival', fontsize=12)
axes[1].set_title('Tumor Suppressor Burden vs PFI', fontsize=11)
axes[1].grid(True, linestyle=':', alpha=0.5)

plt.suptitle('OncoKB Mutation Class Analysis', fontsize=13, fontweight='bold')
plt.tight_layout()
plt.savefig('figures/05_oncokb_mutation_class.png', dpi=300, bbox_inches='tight')
print("\nSaved figures/05_oncokb_mutation_class.png")


# Shoenfeld residuals again for Oncogene
og_cox_data = pfi_data[['PFI.time', 'PFI', 'has_secondary_og']].dropna()
cph_og = CoxPHFitter()
cph_og.fit(og_cox_data, duration_col='PFI.time', event_col='PFI')

ph_test = proportional_hazard_test(cph_og, og_cox_data, time_transform='rank')
print(ph_test.summary)