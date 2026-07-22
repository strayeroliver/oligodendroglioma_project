# This will be mostly copy-pasted from 11, but it would auto-filter in RNA sequences to auto extract the values of expression for the 9 genes.

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from lifelines import CoxPHFitter

oligo = pd.read_csv('data/TCGA/oligo_confirmed_master.csv')
oligo['grade_encoded'] = (oligo['neoplasm_histologic_grade'] == 'G3').astype(int)
oligo['secondary_og_burden'] = oligo[['NOTCH1', 'PIK3CA']].sum(axis=1)
oligo['has_secondary_og'] = (oligo['secondary_og_burden'] > 0).astype(int)

rnaseq = pd.read_csv('data/TCGA/rnaseq_survival_ready.csv')

LASSO_COEFS = {
    'VEPH1': 0.18971, 'TRH': 0.20784, 'HLA-DQA2': 0.11264, 
    'SIX1': 0.04329, 'CENPV': -0.26925, 'ABCC3': 0.25340, 
    'DLX6': 0.22175, 'PAX5': 0.36958, 'SEL1L3': 0.11913
}

rnaseq['rnaseq_risk_score'] = 0.0
for gene, beta in LASSO_COEFS.items():
    if gene in rnaseq.columns:
        rnaseq['rnaseq_risk_score'] += rnaseq[gene] * beta

combined = oligo.merge(rnaseq[['sampleID', 'rnaseq_risk_score']], 
                        on='sampleID', how='inner')
combined = combined[combined['PFI.time'] > 0].dropna(
    subset=['PFI', 'PFI.time', 'grade_encoded', 
            'age_at_initial_pathologic_diagnosis',
            'has_secondary_og', 'rnaseq_risk_score']).copy()
combined['PFI.time'] = combined['PFI.time'] / 365

model_vars = ['PFI.time', 'PFI', 'grade_encoded', 
              'age_at_initial_pathologic_diagnosis',
              'has_secondary_og', 'rnaseq_risk_score']
cox_data = combined[model_vars].dropna()

cph = CoxPHFitter(penalizer=0.1)
cph.fit(cox_data, duration_col='PFI.time', event_col='PFI')

combined['risk_score'] = cph.predict_partial_hazard(
    cox_data[['grade_encoded', 'age_at_initial_pathologic_diagnosis',
              'has_secondary_og', 'rnaseq_risk_score']]
)
training_scores = combined['risk_score'].values

gene_expression_stats = {}
for gene in LASSO_COEFS.keys():
    if gene in rnaseq.columns:
        gene_expression_stats[gene] = {
            'mean': rnaseq[gene].mean(), 'std': rnaseq[gene].std(),
            'min': rnaseq[gene].min(), 'max': rnaseq[gene].max()
        }

def extract_genes_from_file(filepath, sample_column=None):
    # This will accept a full RNA-seq expression file with the same format as the TCGA I downloaded
    # so: genes as rows, patient IDs as columns, all tab seperated and it will auto pull the 9 for the model.
    print(f"Expression file: {filepath}")
    file_data = pd.read_csv(filepath, sep='\t', index_col=0)
    
    if sample_column is None:
        sample_column = file_data.columns[0]
        print(f"No sample specified, using first column: {sample_column}")
    
    patient_series = file_data[sample_column]
    
    extracted = {}
    missing = []
    for gene in LASSO_COEFS.keys():
        if gene in patient_series.index:
            extracted[gene] = patient_series[gene]
        else:
            missing.append(gene)
    
    print(f"{len(extracted)}/9 required genes")
    if missing:
        print(f"Missing genes (excluded from score): {missing}")
    
    return extracted

# And this is all mostly the same

def calculate_rnaseq_score(expression_dict):
    score = 0.0
    contributions = {}
    for gene, beta in LASSO_COEFS.items():
        if gene in expression_dict:
            score += expression_dict[gene] * beta
            cohort_mean = gene_expression_stats.get(gene, {}).get('mean', 0)
            centered_contribution = (expression_dict[gene] - cohort_mean) * beta # This is to actually see how much it deviates from the average
            contributions[gene] = centered_contribution # Previously, i just had "positive coeff = increases risk, neg = lower" which is not true at all, it depends on how diffferent it is compared to the average not to 0.

    return score, contributions
# No need for an else as we already did that earlier

def plot_patient_survival(surv_func, patient_label="Patient"):
    plt.figure(figsize=(10,6))
    plt.plot(surv_func.index, surv_func.values, 
             label=f'Survival Curve ({patient_label})', 
             color='#1f77b4', linewidth=2.5, drawstyle='steps-post')
    plt.plot(cph.baseline_survival_.index, cph.baseline_survival_.values, 
             label='Cohort Baseline Avg', 
             color='gray', linewidth=1.5, linestyle='--', drawstyle='steps-post')
    plt.title(f"Progression-Free Survival Probability for {patient_label}", 
              fontsize=14, fontweight='bold', pad=15)
    plt.xlabel("Time (Years)", fontsize=12)
    plt.ylabel("Survival Probability", fontsize=12)
    plt.ylim(0, 1.05)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend(fontsize=11)
    plt.tight_layout()
    underscored_label = patient_label.replace(" ", "_")
    plt.savefig(f'figures/patients/{underscored_label}_survival_curve.png', dpi=300)
    plt.close() # exact same

def predict_patient(grade_3, age, has_secondary_oncogene, 
                     gene_expression, patient_label="Patient"):
    rnaseq_score, contributions = calculate_rnaseq_score(gene_expression)
    
    patient_df = pd.DataFrame({
        'grade_encoded': [grade_3],
        'age_at_initial_pathologic_diagnosis': [age],
        'has_secondary_og': [has_secondary_oncogene],
        'rnaseq_risk_score': [rnaseq_score]
    }) # same
    
    patient_risk = cph.predict_partial_hazard(patient_df).values[0]
    percentile = (training_scores < patient_risk).mean() * 100
    median_risk = np.median(training_scores)
    risk_class = "HIGH RISK" if patient_risk > median_risk else "LOW RISK"
    surv_func = cph.predict_survival_function(patient_df)
    
    print(f"\nPrediction report for '{patient_label}'")
    print(f"Grade {'III' if grade_3 else 'II'}, Age {age}")
    print(f"Secondary oncogene (NOTCH1/PIK3CA): {'Present' if has_secondary_oncogene else 'Absent'}")
    print(f"\nGene expression contributions to risk score:")
    for gene, contrib in sorted(contributions.items(), key=lambda x: abs(x[1]), reverse=True):
        direction = "increases" if contrib > 0 else "decreases"
        stats = gene_expression_stats.get(gene, {})
        patient_val = gene_expression.get(gene, None)
        cohort_mean = stats.get('mean', None)
        comparison = f"(cohort avg: {cohort_mean:.2f})" if cohort_mean is not None else ""
        print(f"   {gene:10} = {patient_val:.2f} {comparison} --> {direction} risk ({contrib:+.3f})") # Just same as before
    
    print(f"\n transcriptomic score: {rnaseq_score:.3f}")
    print(f"Overall risk classification: {risk_class}")
    print(f"Percentile vs. training cohort (n={len(training_scores)}): {percentile:.1f}th percentile")
    print(f"\nPredicted Progression-Free Survival Probability:")
    for year in [1, 3, 5, 10, 12]:
        prob = surv_func.asof(year).values[0]
        print(f"  {year:2d} years: {prob*100:.1f}%")
    
    plot_patient_survival(surv_func, patient_label)
    return patient_risk, percentile, risk_class, contributions



# If you had a patients data, you would do:
# gene_expr = extract_genes_from_file('data/patient_samples/*******_rnaseq.txt')
# predict_patient(grade_3=#, age=##, has_secondary_oncogene=#, gene_expression=gene_expr, patient_label="*******")




