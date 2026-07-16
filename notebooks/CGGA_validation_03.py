# To validate if what I found in TCGA holds up in CGGA (Chinese cohort of data), or if its just a bunch of baloney!!
# Unfortunately, this dataset ony has OS, not PFI, so the results will be slightly different but we can take it with a grain of salt- this is just for validation, and if the
# trends continue then it would be said to be validated. But, remember, the numbers will not be the exact same.

# https://www.cgga.org.cn/download.jsp
# Download it; https://markgalassi.codeberg.page/small-courses-html/web-scraping/web-scraping.html

# There are two datasets that I will be analyzing inside the CGGA cohort. 325 and 693. 

# Same setup as all the others.

# This python file is essentially the exact same thing as data_loading_01 and cox_regression_02, just shorter, and practically copy-pasted with swapped variables

import matplotlib
matplotlib.use('Agg')  
import pandas as pd
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter, CoxPHFitter

cgga_325 = pd.read_csv('data/CGGA/CGGA.mRNAseq_325_clinical.20200506.txt', sep = '\t')
cgga_693 = pd.read_csv('data/CGGA/CGGA.mRNAseq_693_clinical.20200506.txt', sep = '\t')

# These datasets use CGGA_ID	PRS_type	Histology	Grade	Gender	Age	OS	Censor (alive=0; dead=1)	Radio_status (treated=1;un-treated=0)	Chemo_status (TMZ treated=1;un-treated=0)	IDH_mutation_status	1p19q_codeletion_status	MGMTp_methylation_status
# as their labels, which is rather clunky and non-uniform. https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.rename.html

def clean_cgga_file(df, batch_name):
    df = df.copy()
    df['batch'] = batch_name # This is to preserve independence for merging
    df = df.rename(columns={ # Renaming dataset https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.rename.html
    'Censor (alive=0; dead=1)': 'dead',
    'Radio_status (treated=1;un-treated=0)': 'radiation',
    'Chemo_status (TMZ treated=1;un-treated=0)': 'chemo',
    '1p19q_codeletion_status': 'codel_status',
    'IDH_mutation_status': 'IDH_status',
    'IDH_mut_status': 'IDH_status', # Error from Genome part, essentially its differnet in the genome sequence for litterally no reason...
    'MGMTp_methylation_status': 'MGMT_status'
    })
    return df

cgga_325_renamed = clean_cgga_file(cgga_325, '325')
cgga_693_renamed = clean_cgga_file(cgga_693, '693')
# https://pandas.pydata.org/docs/reference/api/pandas.concat.html to stack them ontop of each other
cgga_merged = pd.concat([cgga_325_renamed, cgga_693_renamed], ignore_index=True) # To count forward from https://pandas.pydata.org/docs/reference/api/pandas.concat.html under 'examples'

print(f"# Of patients in cgga {len(cgga_merged)}") #1018 total!!! success, but this is all glioma- we want oligodendroglioma

cgga_oligo = cgga_merged[
    (cgga_merged['codel_status'] == 'Codel') & # this is isolating the ones with ip19q_codeletion_status, which is essentially needed for oligodendroglioma
    (cgga_merged['IDH_status'] == 'Mutant')  # Again, isolating the one with IDH which is also needed (not making that same mistake again...)
].copy()

print(f"Oligodendroglioma CGGA: {len(cgga_oligo)}") # only 182 oligodendroglioma confirmed
print(f"\nDeaths: {int(cgga_oligo['dead'].sum())}") # and only 50 deaths, which is a pretty small dataset

#Now, we cox it up! I want to track same as before: grade, age, gender, radiation
# We also encode age and mehtylation status to integers, like we did in the last code
cgga_oligo['grade_3plus'] = (cgga_oligo['Grade'].isin(['WHO III', 'WHO IV'])).astype(int)
cgga_oligo['MGMT_methylated'] = (cgga_oligo['MGMT_status'] == 'methylated').astype(int)

cgga_variables = {
    'Grade 3+': 'grade_3plus',
    'Age': 'Age',
    'MGMT Methylated': 'MGMT_methylated',
    'Radiation': 'radiation',
    'Chemotherapy (TMZ)': 'chemo'
} #MGMT is a more-indepth approach for the chemotherapy, and lets you know if the genes are being repaired or not following chemo 

results_cgga=[] #results
cph = CoxPHFitter()
for name, col in cgga_variables.items(): #Univariable approach, just copy pasted from python file 02 essentially.
    try:
        data = cgga_oligo[['OS', 'dead', col]].dropna() # again, drop
        if data[col].nunique() < 2: #same thing for variation
            continue
        cph.fit(data, duration_col='OS', event_col='dead')
        summary = cph.summary
        hr = summary['exp(coef)'].values[0]
        ci_low = summary['exp(coef) lower 95%'].values[0]
        ci_high = summary['exp(coef) upper 95%'].values[0] #confidence intervals
        p_val = summary['p'].values[0]

        results_cgga.append({
            'Variable': name,
            'Hazard Ratio (HR)': round(hr, 3),
            '95% CI': f"{round(ci_low, 2)}-{round(ci_high, 2)}", 
            'p-value': round(p_val, 4)
        })
    except Exception as e:
        print(f"Skipped {name}: {e}")
# This whole thing was transcripted from my journal #2, copied very smoothly. Now, we will copy and paste it again for the biological signatures in CGGA
results_cgga_df = pd.DataFrame(results_cgga).sort_values('p-value')
print("\nCCGA Results for Clinical Values")
print(results_cgga_df.to_string(index=False))
# We got some great results, we see that grade has a extreme HR factor (which is expected) with age following behind. This is as we saw in the previous TCGA file.

# Loading in WEseq dataset: https://www.cgga.org.cn/download.jsp to look at NOTCH1 and other important biological signatures (validation for TCGA, again)
wes_clinical = pd.read_csv('data/CGGA/CGGA.WEseq_286_clinical.20200506.txt', sep='\t')
wes_mutations = pd.read_csv('data/CGGA/CGGA.WEseq_286.20200506.txt', sep='\t') # exact same process

wes_mut_t = wes_mutations.set_index('Unnamed: 0').T # Snippit, https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.transpose.html
wes_mut_t.index.name = 'CGGA_ID' # Snippit
wes_mut_t = wes_mut_t.reset_index() # Gemini snippit. Essentially, it cannot transpose easily as Unnamed:0 is not an index, so the very first column isn't an index meaning it cannot be transposed in the original way unless we make it an index. It will search row numbers, not gene names.
wes_clinical_clean = clean_cgga_file(wes_clinical, 'WES') # Changing variable names like earlier
wes_oligo = wes_clinical_clean[
    (wes_clinical_clean['IDH_status'] == 'Mutant') &
    (wes_clinical_clean['codel_status'] == 'Codel')
].copy() #Filtering for IDH mutated and 1p19q codeletion, again both meaning oligodendroglioma

wes_combined = wes_oligo.merge(wes_mut_t, on='CGGA_ID', how='inner').dropna(subset=['OS', 'dead']) # merging them together, https://www.tutorialspoint.com/article/merge-python-pandas-dataframe-with-a-common-column-and-set-nan-for-unmatched-values
wes_combined['grade_3plus'] = (wes_combined['Grade'].isin(['WHO III', 'WHO IV'])).astype(int) # scanning both grade 3 and 4 (CGGA has grade 4 as well, we are treating them as the same)
wes_combined['male'] = (wes_combined['Gender'] == 'Male').astype(int) # Exact same as we did before
wes_combined['NOTCH1'] = (wes_combined['NOTCH1'].notna()).astype(int)
wes_combined['CIC'] = (wes_combined['CIC'].notna()).astype(int)
wes_combined['FUBP1'] = (wes_combined['FUBP1'].notna()).astype(int)
wes_combined['PIK3CA'] = (wes_combined['PIK3CA'].notna()).astype(int) # Copy and pasted from before

genes_to_test = { # Same COX regression format
    'NOTCH1': 'NOTCH1',
    'CIC': 'CIC',
    'FUBP1': 'FUBP1',
    'PIK3CA': 'PIK3CA'
}

results_cgga_genome = [] 
for name, col in genes_to_test.items(): #Completely copy-pasted from above
    try:
        data = wes_combined[['OS', 'dead', col]].dropna() # again, drop
        if data[col].nunique() < 2: #same thing for variation
            continue
        cph.fit(data, duration_col='OS', event_col='dead')
        summary = cph.summary
        hr = summary['exp(coef)'].values[0]
        ci_low = summary['exp(coef) lower 95%'].values[0]
        ci_high = summary['exp(coef) upper 95%'].values[0] #confidence intervals
        p_val = summary['p'].values[0]

        results_cgga_genome.append({
            'Variable': name,
            'Hazard Ratio (HR)': round(hr, 3),
            '95% CI': f"{round(ci_low, 2)}-{round(ci_high, 2)}", 
            'p-value': round(p_val, 4)
        })
    except Exception as e:
        print(f"Skipped {name}: {e}")
results_cgga_genome_df = pd.DataFrame(results_cgga_genome).sort_values('p-value')
print("\nCCGA Results for Genome Values")
print(results_cgga_genome_df.to_string(index=False))
# The results were completely underpowered. Total of around 50 deaths in the entire set, likely under 10 patients had NOTCH1 and died, very few had PIK3CA and FUBP1. This dataset was way too underpowered for genomic sequence analysis.
# This is because the # of deaths were limited, and hte fact that I had so many biological signatures that were all binary integers being tested really ruined it. The clinical approach was validated, however.
# I can look into more datasets, or I can continue and say TCGA was accurate (and say CGGA was simply underpowered), and continue on my next journey which is part of my larger question
# Are there individual changes in the RNA sequence that promote survival odds? Not just large, biological signatures but the RNA sequence itself.




