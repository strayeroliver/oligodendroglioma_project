# Bootstrap-corrected concordance for the full model.
# This is because my previous concordance of 0.841 might just be overfit. I just tested it on the same 125 patients and 35 events 
# so ofcourse I should get a pretty high probability result. 
# Instead, what we will do is we will fit a model on the bootstrap sample (a resample of the data)
# Then we will score the model on it was on, which will be the same 'apparent concordance' that I found earlier (or same appraoch/mistake)
# Then, score the model on teh original data and just repeat. 
# We will do (concordance bootstrap) - (concordance original) repeatedly, averaging it out, which will give us an apparent understanding of what the overfit concordance is.
# Then just subtract, concordance original (0.841) - average 

import pandas as pd
import numpy as np
from lifelines import CoxPHFitter
from lifelines.exceptions import ConvergenceError
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt




np.random.seed(42)

data = pd.read_csv('data/TCGA/oligo_final_combined_model.csv')
data = data.dropna(subset=['PFI.time', 'PFI'] + ['grade_encoded', 'age_at_initial_pathologic_diagnosis',
              'has_secondary_og', 'rnaseq_risk_score', 'PAX5_methylation']).copy()
model_cols = ['grade_encoded', 'age_at_initial_pathologic_diagnosis',
              'has_secondary_og', 'rnaseq_risk_score', 'PAX5_methylation'] + ['PFI.time', 'PFI']



def fit_cox(df): # I wanted to try out a new function style to see if it is easier for typing, and just to learn more about functions and lifelines
    cph = CoxPHFitter(penalizer=0.1)
    cph.fit(df[model_cols], duration_col='PFI.time', event_col='PFI')
    return cph


cph_full = fit_cox(data)
apparent_c = cph_full.score(data[model_cols], scoring_method='concordance_index')
print(f"\nApparent concordance: {apparent_c:.3f}")



# Now onto bootstraps
optimism_estimates = []
corrected_c_estimates = []
n_failed = 0

for b in range(500): # 500 bootstraps more the merrier
    boot_sample = data.sample(n=len(data), replace=True, random_state=42 + b)

    # skip resamples that don't have enough variation/events to converge 
    # this happens sometimes with small cohorts and is normal, just log it
    if boot_sample['PFI'].sum() < 5:
        n_failed += 1
        continue

    try:
        cph_boot = fit_cox(boot_sample)
    except ConvergenceError:
        n_failed += 1
        continue

    c_boot_train = cph_boot.score(boot_sample[model_cols], scoring_method='concordance_index')
    c_boot_on_orig = cph_boot.score(data[model_cols], scoring_method='concordance_index')

    optimism_b = c_boot_train - c_boot_on_orig
    optimism_estimates.append(optimism_b)
    corrected_c_estimates.append(apparent_c - optimism_b)

if n_failed > 0:
    print(f"({n_failed}/{500} bootstrap resamples skipped due to convergence/lack of events)")

optimism_estimates = np.array(optimism_estimates)
corrected_c_estimates = np.array(corrected_c_estimates)

mean_optimism = optimism_estimates.mean()
corrected_c = apparent_c - mean_optimism # as i said earlier


ci_low, ci_high = np.percentile(corrected_c_estimates, [2.5, 97.5])

print(f"\nBootstrap (mean over {len(optimism_estimates)}): {mean_optimism:.3f}")
print(f"Optimism-corrected concordance: {corrected_c:.3f}")
print(f"95% CI (bootstrap percentile): [{ci_low:.3f}, {ci_high:.3f}]")

plt.figure(figsize=(8, 5))
plt.hist(corrected_c_estimates, bins=30, color='#2c7bb6', edgecolor='black', alpha=0.8)
plt.axvline(apparent_c, color='#d7191c', linestyle='--', linewidth=2,
            label=f'Apparent c = {apparent_c:.3f}')
plt.axvline(corrected_c, color='orange', linestyle='-', linewidth=3,
            label=f'Corrected c = {corrected_c:.3f}')
plt.xlabel('Bootstrap-corrected concordance')
plt.ylabel('Frequency')
plt.title(f'Bootstrap Optimism Correction (n={len(corrected_c_estimates)} valid resamples)',
          fontsize=13, fontweight='bold')
plt.legend()
plt.grid(True, linestyle=':', alpha=0.5)
plt.tight_layout()
plt.savefig('figures/14_bootstrap_correction.png', dpi=300, bbox_inches='tight')
