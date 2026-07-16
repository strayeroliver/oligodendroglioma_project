# This will be the same as rna_analysis but will be using CCGA for validation, hopefully replicating the results. 
# LASSO will NOT NOT NOT be need to be done again, we just transfer over the variables
# Thus, a really quick notebook, should be very fast as it is the exact same as previous notebooks with different variables.

import matplotlib
matplotlib.use('Agg')  
import pandas as pd
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines import CoxPHFitter

# Copy-pasted from earlier
cgga_325_clin = pd.read_csv('data/CGGA/CGGA.mRNAseq_325_clinical.20200506.txt', sep='\t')
cgga_693_clin = pd.read_csv('data/CGGA/CGGA.mRNAseq_693_clinical.20200506.txt', sep='\t')
# Load in RNA seq data also downloaded from same website
cgga_rna_325 = pd.read_csv('data/CGGA/CGGA.mRNAseq_325.RSEM-genes.20200506.txt', sep='\t', index_col=0)
cgga_rna_693 = pd.read_csv('data/CGGA/CGGA.mRNAseq_693.RSEM-genes.20200506.txt', sep='\t', index_col=0)

# Merge
cgga_clin_merged = pd.concat([cgga_325_clin, cgga_693_clin], ignore_index=True)

# search for 1p19 codel and idh mutated for oligodendroglioma molecular confirmed 
oligo_patients = cgga_clin_merged[
    (cgga_clin_merged['1p19q_codeletion_status'] == 'Codel') & 
    (cgga_clin_merged['IDH_mutation_status'] == 'Mutant')
].copy()
oligo_ids = oligo_patients['CGGA_ID'].tolist() # isolate the patient IDs for future use

# Done before in different file. Just filtering out purely from the oligo_ids file to get purely rna sequences of oligodendroglioma patients
cgga_rna_325_filtered = cgga_rna_325[[id for id in oligo_ids if id in cgga_rna_325.columns]]
cgga_rna_693_filtered = cgga_rna_693[[id for id in oligo_ids if id in cgga_rna_693.columns]]

# Merge RNA and transpose. Same process as other codes.
cgga_rna_combined = pd.concat([cgga_rna_325_filtered, cgga_rna_693_filtered], axis=1)
cgga_rna_transposed = cgga_rna_combined.T
cgga_rna_transposed.index.name = 'CGGA_ID'
cgga_rna_transposed = cgga_rna_transposed.reset_index()

# Rename the result info and merge them together
oligo_outcomes = oligo_patients[['CGGA_ID', 'OS', 'Censor (alive=0; dead=1)']].rename(
    columns={'Censor (alive=0; dead=1)': 'dead'}
)
cgga_rna_df = cgga_rna_transposed.merge(oligo_outcomes, on='CGGA_ID', how='inner')
cgga_rna_df = cgga_rna_df.dropna(subset=['OS', 'dead'])
cgga_rna_df = cgga_rna_df[cgga_rna_df['OS'] > 0] # All this is transcripted directly just renamed variables

# from rna_analysis_04.py and R file

coefs = {
    'VEPH1': 0.18971, 'ADAM6': 0.11343, 'TRH': 0.20784, 'HLA-DQA2': 0.11264, 
    'SIX1': 0.04329, 'CENPV': -0.26925, 'ABCC3': 0.25340, 'DLX6': 0.22175, 
    'C18orf34': -0.27209, 'PAX5': 0.36958, 'SEL1L3': 0.11913
}

# Same as KMF stuff
cgga_rna_df['Score'] = 0.0 
for gene, beta in coefs.items():
    if gene in cgga_rna_df.columns:
        cgga_rna_df['Score'] += cgga_rna_df[gene] * beta
    else:
        print(f"{gene} missign") # We see that 2 genes are missing, ADAM6 and C18orf34. We can move on for now as they are fully not inside the dataset after looking inside the actual file.

median_value = cgga_rna_df['Score'].median()
cgga_rna_df['HighRisk'] = cgga_rna_df['Score'] > median_value
cgga_rna_df['LowRisk'] = cgga_rna_df['Score'] <= median_value


fig, ax = plt.subplots(figsize=(10, 6))
kmf = KaplanMeierFitter()

groups = [
    (cgga_rna_df['HighRisk'], 'High Risk Group', '#E41A1C'),
    (cgga_rna_df['LowRisk'], 'Low Risk Group', '#377EB8')
]

for mask, name, color in groups:
    kmf.fit(cgga_rna_df.loc[mask, 'OS'], cgga_rna_df.loc[mask, 'dead'], label=name)
    kmf.plot_survival_function(ax=ax, color=color, lw=2.5, ci_show=True) 

plt.title('OS CGGA RNA Seperation', pad=12, fontsize=14, weight='bold')
plt.xlabel('OS Time in Days', fontsize=12)
plt.ylabel('Survival Probability', fontsize=12)
plt.grid(True, linestyle=':', alpha=0.5)

output_path = 'figures/KM_CCGA_validation.png'
plt.savefig(output_path, dpi=300, bbox_inches='tight')


# To check if they are actually seperated via genes instead of just "Oh we took a random bunch of people, split them in the median, and hey look they are split!"
# We have to see Log-Rank Test: https://biostatsquid.com/easy-log-rank-test/
# https://medium.com/@leushina.katya/nonparametric-methods-in-survival-analysis-log-rank-test-c1bda672d7b9
# https://pmc.ncbi.nlm.nih.gov/articles/PMC403858/

from lifelines.statistics import logrank_test 
# Careful help from https://medium.com/@leushina.katya/nonparametric-methods-in-survival-analysis-log-rank-test-c1bda672d7b9
results = logrank_test(
    cgga_rna_df.loc[cgga_rna_df['HighRisk'], 'OS'], cgga_rna_df.loc[cgga_rna_df['LowRisk'], 'OS'],
    event_observed_A=cgga_rna_df.loc[cgga_rna_df['HighRisk'], 'dead'], event_observed_B=cgga_rna_df.loc[cgga_rna_df['LowRisk'], 'dead']
)

print(f"Log-Rank Test p-value: {results.p_value:.6f}")
# P value of 0.029677 < 0.05, so this means this is not just cutting in half randomly it means there is actually significant difference between the two groups.

# To find concordance of model so we rerun cox, not just log-rank test for high/low risk.
# Same idea as what we did before
gene_cols = [g for g in coefs.keys() if g in cgga_rna_df.columns] # Extract columns of genes
cox_data_cgga = cgga_rna_df[['OS', 'dead'] + gene_cols].copy() # Attach the genes to the clinical data (death etc)
cox_data_cgga = cox_data_cgga[cox_data_cgga['OS'] > 0] # Make sure all OS>0
cox_data_cgga = cox_data_cgga.dropna() # Drop missing values, exact same as we have done in the past.
cox_data_cgga[gene_cols] = (cox_data_cgga[gene_cols] - cox_data_cgga[gene_cols].mean()) / cox_data_cgga[gene_cols].std() # https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.std.html, new standardized values
cox_cgga = CoxPHFitter(penalizer=0.1) # This time we attach a penalizer because we are testing p > n (a very wishy-washy lasso-esque idea)
cox_cgga.fit(cox_data_cgga, duration_col='OS', event_col='dead')
cox_cgga.print_summary(decimals=3) # We get a concordance of 0.734 which is good so we are dancing in the streets!!

