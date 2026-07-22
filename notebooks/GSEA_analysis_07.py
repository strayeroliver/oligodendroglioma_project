# Gene Set Enrichment Analysis (GSEA) - pathway-level interpretation of the RNA-seq risk score findings

# This looks at why the LASSO genes that we found earlier in file 04 actually mean what they mean. 
# I briefly looked at some of them online and saw a few associated with T-cells/immune system, all things surrounding the brain, etc.
# Now we look at their 'classification' in a pathway. For example, it would classify the genes into pathway-specific labels.
# It takes thousands of genes and asks which biological programs are running differently between high and low risk patients
# trying to see the larger pathway thats changing due to high vs low risk patients.

# If a large proportion of genes in a pathway are changed and it shows increase in risk, then that pathway probably links to something.

# I will use gseapy (Python implementation of Broad Institute GSEA)
# https://gseapy.readthedocs.io/en/latest/gseapy_example.html for GSEAPY help

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy import stats
from lifelines import CoxPHFitter
import gseapy as gp
import os

oligo = pd.read_csv('data/TCGA/oligo_confirmed_master.csv')
rnaseq_survival = pd.read_csv('data/TCGA/rnaseq_survival_ready.csv')
rnaseq_survival = rnaseq_survival[rnaseq_survival['PFI.time'] > 0]

# LASSO selected genes and coefficients from R output
# These are the mean-variance filtered LASSO results (11 genes)
lasso_coefs = {
    'VEPH1': 0.18971, 'ADAM6': 0.11343, 'TRH': 0.20784,
    'HLA-DQA2': 0.11264, 'SIX1': 0.04329, 'CENPV': -0.26925,
    'ABCC3': 0.25340, 'DLX6': 0.22175, 'C18orf34': -0.27209,
    'PAX5': 0.36958, 'SEL1L3': 0.11913
}

# Copy pasted from previous file
rnaseq_survival['risk_score'] = 0.0
for gene, beta in lasso_coefs.items():
    if gene in rnaseq_survival.columns:
        rnaseq_survival['risk_score'] += rnaseq_survival[gene] * beta

median_score = rnaseq_survival['risk_score'].median()
rnaseq_survival['risk_group'] = (rnaseq_survival['risk_score'] > median_score).astype(int)


# Full 5000-gene filtered RNA-seq matrix done in rna_analysis_04 - to actually see the larger amount for pathways 

rnaseq_full = pd.read_csv('data/TCGA/rnaseq_survival_ready.csv')
rnaseq_full = rnaseq_full[rnaseq_full['PFI.time'] > 0]

# Merging both full and survival to get risk groups 
rnaseq_full = rnaseq_full.merge(
    rnaseq_survival[['sampleID', 'risk_group']],
    on='sampleID', how='inner'
)
# setting gene columns 
gene_cols = [c for c in rnaseq_full.columns
             if c not in ['sampleID', 'PFI', 'PFI.time', 'risk_group']]

# Rank by t-statistic: high risk vs low risk # https://en.wikipedia.org/wiki/T-statistic
# Positive t-stat = higher expression in high-risk patients # https://blog.minitab.com/en/blog/statistics-and-quality-data-analysis/what-are-t-values-and-p-values-in-statistics
# Negative t-stat = higher expression in low-risk patients
high_risk = rnaseq_full[rnaseq_full['risk_group'] == 1]
low_risk = rnaseq_full[rnaseq_full['risk_group'] == 0]

ranking = {}
for gene in gene_cols:
    h_vals = high_risk[gene].dropna()
    l_vals = low_risk[gene].dropna()
    if len(h_vals) > 5 and len(l_vals) > 5:
        t_stat, _ = stats.ttest_ind(h_vals, l_vals) # Using equation to find t-statistic
        ranking[gene] = t_stat

ranked_genes = pd.Series(ranking).sort_values(ascending=False) # https://pandas.pydata.org/docs/reference/api/pandas.Series.html makes it a column and sorts by size
print(f"Ranked genes: {len(ranked_genes)}")
print(f"\nTop 5 high risk:")
print(ranked_genes.head())
print(f"\nTop 5 low risk:")
print(ranked_genes.tail())


# GSEA has three databases, we want to essentially test our 5000 genes on these:
# https://gseapy.readthedocs.io/en/latest/gseapy_example.html for all gp. commands and learning about it

os.makedirs('results/gsea_hallmark', exist_ok=True) # hallmark, kegg, reactome are the three databases
os.makedirs('results/gsea_kegg', exist_ok=True)
os.makedirs('results/gsea_reactome', exist_ok=True)

# https://www.genepattern.org/modules/docs/GSEAPreranked/1/
print("\nHallmark")
hallmark = gp.prerank(
    rnk=ranked_genes,
    gene_sets='MSigDB_Hallmark_2020',
    threads=4,
    min_size=15, # this is min size of pathway, so pathways less than 15 genes and larger than 500 are excluded.
    max_size=500,
    permutation_num=1000, # default
    outdir='results/gsea_hallmark',
    seed=42,
    verbose=False
)
hallmark_df = hallmark.res2d.sort_values('NES', ascending=False) 
sig_hallmark = hallmark_df[hallmark_df['FDR q-val'] < 0.25] # essentially its a p value, anything less than 0.25 is reasonably sound
print(f"Significant Hallmark pathways: {len(sig_hallmark)}")  


print("KEGG")
kegg = gp.prerank(
    rnk=ranked_genes,
    gene_sets='KEGG_2021_Human',
    threads=4,
    min_size=15,
    max_size=500,
    permutation_num=1000,
    outdir='results/gsea_kegg',
    seed=42,
    verbose=False # Exact same as the prior
)
kegg_df = kegg.res2d.sort_values('NES', ascending=False)
sig_kegg = kegg_df[kegg_df['FDR q-val'] < 0.25]
print(f"Significant KEGG pathways: {len(sig_kegg)}")

print("Reactome")
reactome = gp.prerank(
    rnk=ranked_genes,
    gene_sets='Reactome_2022',
    threads=4,
    min_size=15,
    max_size=500,
    permutation_num=1000,
    outdir='results/gsea_reactome',
    seed=42,
    verbose=False
)
reactome_df = reactome.res2d.sort_values('NES', ascending=False)
sig_reactome = reactome_df[reactome_df['FDR q-val'] < 0.25]
print(f"Significant Reactome pathways: {len(sig_reactome)}")



# Now we summarize and output
# NES scores are essentially how prominent they are or how much they matter (essentially HR from previous examples)

print("\nHALLMARK")
print(sig_hallmark[['Term', 'NES', 'NOM p-val', 'FDR q-val']].head(10).to_string(index=False))  # All sorted by NES

print("\nKEGG (primarily oncology/cancer related dataset)")
cancer_kegg = kegg_df[kegg_df['Term'].str.contains(
    'Cell cycle|PI3K|Notch|MAPK|p53|Cancer|Wnt|mTOR', # filtering
    case=False, na=False
)]
print(cancer_kegg[['Term', 'NES', 'NOM p-val', 'FDR q-val']].head(10).to_string(index=False))
print(sig_kegg[['Term', 'NES', 'NOM p-val', 'FDR q-val']].head(10).to_string(index=False)) # unfiltered as when I ran earlier it said nothing was significant

print("\nREACTOME (primarily cell cycle and division dataset)")
cell_cycle_react = reactome_df[reactome_df['Term'].str.contains(
    'Cell Cycle|Notch|PI3K|Mitotic|Checkpoint', # filtering
    case=False, na=False
)]
print(cell_cycle_react[['Term', 'NES', 'NOM p-val', 'FDR q-val']].head(10).to_string(index=False))

# We see a lot of different pathways and values/results.
# HALLMARK, which is the more common dataset others use, provides the best results + pathways so we will be making a graph with hallmark only


sig_sorted = sig_hallmark.sort_values('NES', ascending=True)

# https://matplotlib.org/stable/gallery/lines_bars_and_markers/barh.html for bar graph help

fig, ax = plt.subplots(figsize=(10, max(8, len(sig_sorted) * 0.3))) # scales based on sig bars
colors = ['#d7191c' if nes > 0 else '#2c7bb6' for nes in sig_sorted['NES']] # the graphs came out all red which makes sense
ax.barh(sig_sorted['Term'], sig_sorted['NES'], color=colors)
ax.axvline(x=0, color='black', linewidth=0.8)
ax.set_xlabel('Normalized Enrichment Score (NES)', fontsize=12)
ax.set_title('GSEA Hallmark Pathways\nHigh Risk vs Low Risk Oligodendroglioma (FDR < 0.25)',
             fontsize=12, fontweight='bold')

# Gemini snippit for how on earth you add labels to a bar graph (like on the bar itself, not on the sides)
for i, (nes, fdr) in enumerate(zip(sig_sorted['NES'], sig_sorted['FDR q-val'])): # I understand now, this is just making it long/short
    ax.text(nes + 0.02 if nes > 0 else nes - 0.02, i, # This is adding text itself about the FDR score 
            f'FDR={fdr:.3f}', va='center', fontsize=7,
            ha='left' if nes > 0 else 'right') 

plt.tight_layout()
plt.savefig('figures/06_gsea_hallmark.png', dpi=300, bbox_inches='tight')
# Just saving results

hallmark_df.to_csv('results/gsea_hallmark_results.csv', index=False)
kegg_df.to_csv('results/gsea_kegg_results.csv', index=False)
reactome_df.to_csv('results/gsea_reactome_results.csv', index=False)

# Creating my full summary of the data
summary = pd.DataFrame({
    'Database': ['Hallmark', 'Hallmark', 'Hallmark', 'Hallmark',
                 'KEGG', 'Reactome', 'Reactome'],
    'Pathway': ['G2-M Checkpoint', 'TNF-alpha Signaling via NF-kB', 
                'Inflammatory Response', 'EMT',
                'Human T-cell Leukemia Virus 1 Infection', 
                'Cell Cycle Mitotic', 'Cell Cycle'],
    'NES': [2.034, 1.829, 1.768, 1.671, 2.095, 2.175, 2.023],
    'FDR': [0.000, 0.010, 0.018, 0.029, 0.003, 0.001, 0.003],
    'Biological Meaning': [ # https://www.gsea-msigdb.org/gsea/msigdb ALL MEANINGS ARE FOUND FROM HERE
        'Cell cycle checkpoint override - cells dividing uncontrolled',
        'Inflammatory/innate immune signaling activated',
        'General inflammatory pathway activation',
        'Loss of neural identity - tumors becoming more invasive',
        'Gene set captures general immune signaling machinery (not literal viral infection)',
        'Mitotic cell cycle activation confirmed in Reactome',
        'Cell cycle activation confirmed in Reactome (second gene set)'
    ] # just labeling for future use also written in notebook
})

summary.to_csv('results/gsea_cross_database_summary.csv', index=False)
print("done")

