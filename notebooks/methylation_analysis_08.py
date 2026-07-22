# This file will be the most computationally intensive and probably where I struggled the most in writing. 
# It will test whether DNA methylation analysis adds predictive value beyond RNA and mutations.
# https://www.illumina.com/techniques/multiomics/epigenetics/dna-methylation-analysis.html
# It is a large 450k array 

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import pandas as pd
import numpy as np
import json
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
from scipy import stats

oligo = pd.read_csv('data/TCGA/oligo_confirmed_master.csv')
confirmed_ids = oligo['sampleID'].tolist()

# Methylation Data loading (Also from UCSC Xena data website, its 450k and gene mapping + methylation data)
methylation = pd.read_csv('data/TCGA/LGG_methylation450k.gz',
                           sep='\t', compression='gzip', index_col=0)
print(f"Full methylation: {methylation.shape}") # Checking to make sure it is the correct one I downloaded



# Load probe-gene mapping, this is from the second link under methylation data (it tells us what genes specifically map to each methylation data as they dont go by same gene names)
# Right now, this file is a huge mess and needs to be written in something in common with all my datasets: the gene names specifically
# Right now, the names are in format 'cg#######' which is not good. Instead, we need to make the first column be genes and have the rest of them be like an 
# encyclopedia with 'cg#####', 'cg#######' to be the rest of the columns to let us know which ones they are near. 
# Furthermore, I downloaded https://support.illumina.com/downloads/infinium_humanmethylation450_product_files.html
# The Manifest File to tell me promoter regions/distance from promoter regions to eliminate body noise.
# Now, the probes and methylation data actually matters instead of just being all throughout.

manifest_path = 'data/TCGA/humanmethylation450_15017482_v1-2.csv'

# The file has 7 rows of manufacturing data so we skip 7, https://stackoverflow.com/questions/24251219/pandas-read-csv-low-memory-and-dtype-options after erroring out on with low memory dtype error
manifest = pd.read_csv(manifest_path, sep=',', skiprows=7, low_memory=False)

# theres a lot of white space in columns, continuous errors: https://pandas.pydata.org/docs/reference/api/pandas.Series.str.strip.html
manifest.columns = manifest.columns.str.strip() 

key_genes = ['DLX6', 'SIX1', 'TRH', 'ABCC3', 'PAX5', 'CENPV', 'MGMT', 'CDKN2A']
promoter_probes = {}
print("\nPromoter Site Probes")
for gene in key_genes:
   
    gene_mask = manifest['UCSC_RefGene_Name'].str.contains(rf'\b{gene}\b', na=False, regex=True) # Look for my gene ONLY (not like DLX61, only DLX6)
    
    # Mask B: Isolate regulatory switches (TSS200, TSS1500, 5'UTR, 1stExon)
    promoter_mask = manifest['UCSC_RefGene_Group'].str.contains('TSS200|TSS1500|5\'UTR|1stExon', na=False, regex=True) # TSS is distance from promoter, filtering out those who are further
    
    # Take Probe IDS fulfilling both 
    matching_ids = manifest[gene_mask & promoter_mask]['IlmnID'].tolist()
    
    # and then we take the ones with data in both
    valid_promoter_probes = [p for p in matching_ids if p in methylation.index]
    
    # save
    promoter_probes[gene] = valid_promoter_probes
    print(f"{gene}: {len(valid_promoter_probes)} probes promoter")

# Renamed
gene_probes = promoter_probes

# Save for future use now
with open('data/TCGA/key_gene_probes.json', 'w') as f:
    json.dump(gene_probes, f, indent=2)


with open('data/TCGA/key_gene_probes.json', 'r') as f:
    gene_probes = json.load(f)


# Filter to patietns purely with methylation data
ids_in_meth = [id for id in confirmed_ids if id in methylation.columns]
print(f"methylation patients: {len(ids_in_meth)}") # We have all 127 patients inside here with methylation data
meth_oligo = methylation[ids_in_meth] 

gene_methylation = {}
for gene, probes in gene_probes.items():
    valid_probes = [p for p in probes if p in meth_oligo.index] # searches for probes in the oligo patient ids, if theres no probes for that gene it skips
    if len(valid_probes) == 0:
        print(f"{gene}: no probes in dataset, skipping")
        continue
    gene_avg = meth_oligo.loc[valid_probes].mean(axis=0)
    gene_methylation[gene] = gene_avg # We have to average because genes have a TON of probes all throughout the gene. We need to average them to see a total of what extent they
                                        # methylation level is. Also, overfitting for COX model if we have way too many variables.

# Scale beta values to percentage points (0-100) instead of proportions (0-1)
# A HR on the 0-1 scale represents going from 0% to 100% methylated in one step,  which produces artificially huge hazard ratios
for gene in gene_methylation.keys():
    gene_methylation[gene] = gene_methylation[gene] * 100

meth_genes_df = pd.DataFrame(gene_methylation) # dataframe with genes and just the averages
meth_genes_df.index.name = 'sampleID'
meth_genes_df = meth_genes_df.reset_index() # Preparing to merge


# Merge
oligo_meth = oligo.merge(meth_genes_df, on='sampleID', how='inner')
oligo_meth = oligo_meth[oligo_meth['PFI.time'] > 0].dropna(
    subset=['PFI', 'PFI.time']
).copy() # same as always
oligo_meth['PFI.time'] = oligo_meth['PFI.time'] / 365 # standardizing years

print(f"\nPatients w/ both clinical + meth {len(oligo_meth)}") # Only 125 so we lost 2 but still good


# COX regression on methylation data + PFI 

gene_cols = [g for g in gene_probes.keys() if g in oligo_meth.columns]

cph = CoxPHFitter()
meth_results = []

for gene in gene_cols:
    data = oligo_meth[['PFI.time', 'PFI', gene]].dropna()
    if len(data) < 20: # we want larger sample sizes
        continue
    try:
        cph.fit(data, duration_col='PFI.time', event_col='PFI')
        s = cph.summary
        hr = s['exp(coef)'].values[0]
        p = s['p'].values[0]
        meth_results.append({
            'Gene': gene,
            'HR': round(hr, 3),
            'p-value': round(p, 4),
            'N': len(data) # This is to check how many 
        })
    except Exception as e:
        print(f"{gene}: error - {e}")

meth_results_df = pd.DataFrame(meth_results).sort_values('p-value')
print(meth_results_df.to_string(index=False))
meth_results_df.to_csv('results/methylation_gene_cox.csv', index=False)


# In total, we see that PAX5 stays on top with HR 1.06, p 0.001. This means for every 1% increase in methylation for PAX5, 6% increase hazard. 
# Also, CENPV is on the verge of that too. This means that hypermethylation of these two increases the risk significantly.

# We also see MGMT is 0.7734 P value, which was confusing at first but after some reasearch:
# MGMT is so common, so like the IDH/1p19q misclassification at the cox_regression_02.py, its so common in all of them that cox cannot read the signal correctly.
# Even though IDH was so common in the oligodendrogliomas, there was only a 57% response back.



