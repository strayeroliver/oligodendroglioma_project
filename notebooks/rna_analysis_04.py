# This will be the python file for HOPEFULLY successful full RNA-sequence analysation using cox-regression and LASSO approach. 
# Goal: Hopefully I find specific, small genes that contribute to high risk rather than purely looking at major biological signatures (NOTCH1, etc)


import matplotlib
matplotlib.use('Agg')  
import pandas as pd
import matplotlib.pyplot as plt
from lifelines import CoxPHFitter

# Data from: https://xenabrowser.net/datapages/?cohort=TCGA%20Lower%20Grade%20Glioma%20(LGG)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
# Download gene expression RNAseq IlluminaHiSeq

rnaseq = pd.read_csv('data/TCGA/LGG_RNAseq.gz', # It is a zipped file so we need to unzip, also the first column in the file is the gene itself so we set that to index
                      sep='\t',
                      compression='gzip',
                      index_col=0)
oligo_confirmed = pd.read_csv('data/TCGA/oligo_confirmed_master.csv')

# Time to filter down this massive dataset... like 76 megabytes of just text. However, the RNA-seq dataset is like 20000x500, we only want the ~127 with oligodendroglioma
ids = oligo_confirmed['sampleID'].tolist() # This is extracting the patient IDs to match with RNA-seq data
ids_rnaseq = [id for id in ids if id in rnaseq.columns] # Matching that the id in oligo_confirmed and ids_rnaseq
print(f"{len(ids_rnaseq)}") # SUCCESS we got a grand total of..... 127 oligodendroglioma rna-seq patients, so now we just have to cut the others
rnaseq_oligo = rnaseq[ids_rnaseq]
print(f"{len(rnaseq_oligo)}") # 20530 total genes, a lot of these are just clutter. I want to take the top 5000 variated genes, then use cox-lasso to get to a top 10.





# Research on how to filter using pandas and etc dropping low-variance genes: https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.filter.html
# https://www.biostars.org/p/481351/, currently genes are columns so:
# https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.var.html#pandas.DataFrame.var was most crucial here
# gene_variance = rnaseq_oligo.var(axis=1) # axis for looking at columns 
# print(gene_variance.describe()) # We get a list of 20530 numbers from 0 to 1, so it looks like it was successful. Now we take the most variated 5000, which are the largest 5000 numbers
# variated_genes = gene_variance.nlargest(5000).index # pandas.pydata.org/docs/reference/api/pandas.DataFrame.nlargest.html
# rnaseq_filtered = rnaseq_oligo.loc[variated_genes] # https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.loc.html#pandas.DataFrame.loc
# Now, we have a dataset of left columns being genes, all the top 5000 most variated genes. We need to transpose then merge with PFI/OS time stuff to do cox
# rnaseq_transposed = rnaseq_filtered.T # https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.transpose.html
# rnaseq_transposed.index.name = 'sampleID' # Labeling this as index, then moving it in the next line as atlas
# rnaseq_transposed = rnaseq_transposed.reset_index() # https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.reset_index.html
# rnaseq_survival = rnaseq_transposed.merge(oligo_confirmed[['sampleID', 'PFI', 'PFI.time']], on='sampleID', how='inner') # classic merge
# rnaseq_survival = rnaseq_survival.dropna(subset=['PFI', 'PFI.time']) 
# rnaseq_survival.to_csv('data/TCGA/rnaseq_survival_ready.csv', index=False) # Saving it for R visualization/LASSO regression.

#### This upwards block was commented out because 5000 genes was still too much data for my set. I have the option to truncate it further to 1000, but purely cutting
#### values off based on their variation doesnt make too much sense (it does for larger n, but 20000 to 1000 purely off that seems much)
#### Instead, I will truncate my dataset using 'Mean-Variance', where we first look at genes that are expressed a decent amount, then we look at their variation.
#### This makes it so we dont have random, useless genes that just happen to be variated, but actually do something.
#### https://www.biostars.org/p/480419/ mean variance idea




expressed_genes = rnaseq_oligo[rnaseq_oligo.mean(axis=1) > 2.0] # gemini snippit help, understood with https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.mean.html, avg expression level, >2.0 as my data is already log2(n) from TCGA website
gene_variance = expressed_genes.var(axis=1) # and now we go back to what we did prior
variated_genes = gene_variance.nlargest(1000).index # copy pasted from above but with 1000
rnaseq_filtered = expressed_genes.loc[variated_genes]
print(gene_variance.describe()) # We see that 16014 genes passed the expressed_genes barrier, then they are all varying extremely following the .var command. Looks good

# Following is copy and pasted from above
rnaseq_transposed = rnaseq_filtered.T # https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.transpose.html
rnaseq_transposed.index.name = 'sampleID' # Labeling this as index, then moving it in the next line as atlas
rnaseq_transposed = rnaseq_transposed.reset_index() # https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.reset_index.html

rnaseq_survival = rnaseq_transposed.merge(oligo_confirmed[['sampleID', 'PFI', 'PFI.time']], on='sampleID', how='inner') # classic merge
rnaseq_survival = rnaseq_survival.dropna(subset=['PFI', 'PFI.time'])
rnaseq_survival = rnaseq_survival[rnaseq_survival['PFI.time'] > 0]
rnaseq_survival.to_csv('data/TCGA/rnaseq_survival_ready.csv', index=False) # Saving it for R visualization/LASSO regression.


# We have to use 'LASSO' "Least Absolute Shrinkage and Selection Operator"
# Adds a penalty to the function (lambda) and shrinks coefficients of non-important variables. 
# Makes it significantly easier to read and understand at the cost of slight chance of shrinking a major gene. 
# Doing it in R because it is simpler, and I want to learn R, and 'Victoria' sits at my table and she can help me as she is learning R deeply.
# https://www.serdarbalci.com/jsurvival/articles/09-lassocox-comprehensive.html https://glmnet.stanford.edu/articles/Coxnet.html
# https://www.r-bloggers.com/2020/05/quick-tutorial-on-lasso-regression-with-example/ , https://www.youtube.com/watch?v=5GZ5BHOugBQ

# Just did LASSO/COX in R, got fantastic results! 
# Now, I will visualize in KMF graph

from lifelines import KaplanMeierFitter # Adding this in

rnaseq_survival = rnaseq_survival[rnaseq_survival['PFI.time'] > 0]

coefs = {
    'VEPH1': 0.18971, 'ADAM6': 0.11343, 'TRH': 0.20784, 'HLA-DQA2': 0.11264, 
    'SIX1': 0.04329, 'CENPV': -0.26925, 'ABCC3': 0.25340, 'DLX6': 0.22175, 
    'C18orf34': -0.27209, 'PAX5': 0.36958, 'SEL1L3': 0.11913
} # This is directly from my R output. The R coefficients, will be used to calculate risk score to then be plotted

rnaseq_survival['Score'] = 0.0
for gene, beta in coefs.items():
    rnaseq_survival['Score'] += rnaseq_survival[gene] * beta # summing up all the scores to get a total
median_value = rnaseq_survival['Score'].median()
rnaseq_survival['HighRisk'] = rnaseq_survival['Score'] > median_value # This is to split accross the median into a 50/50 high/low risk
rnaseq_survival['LowRisk'] = rnaseq_survival['Score'] <= median_value

fig, ax = plt.subplots(figsize=(10, 6))
kmf = KaplanMeierFitter()

groups = [
    (rnaseq_survival['HighRisk'], 'High Risk Group', '#E41A1C'),
    (rnaseq_survival['LowRisk'], 'Low Risk Group', '#377EB8')
]

for mask, name, color in groups:
    kmf.fit(rnaseq_survival.loc[mask, 'PFI.time'], rnaseq_survival.loc[mask, 'PFI'], label=name)
    kmf.plot_survival_function(ax=ax, color=color, lw=2.5, ci_show=True) # This for loop was done by gemini. I understand it fully, now. Essentially, (below)
    # It takes in the data (high/low, name, red/blue) and cuts in half using .loc[mask]. It fits the graph, then plots it under.

plt.title('Oligodendroglioma Progression-Free Interval', pad=12, fontsize=14, weight='bold')
plt.xlabel('Time in Days', fontsize=12)
plt.ylabel('Progression-Free Probability', fontsize=12)
plt.grid(True, linestyle=':', alpha=0.5)

output_path = 'figures/km_LASSO.png'
plt.savefig(output_path, dpi=300, bbox_inches='tight')





































