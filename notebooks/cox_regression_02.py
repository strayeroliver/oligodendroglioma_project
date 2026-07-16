# For actual clinical variables now- age, gender, etc, and also major biological mutations from previous research done on Oligodendroglioma

# Although this python script looks long, it is really the same thing just repeated 4 times. I just made a mistake with the dataset and results and had to clean it, but left the old cox-regressions in there.
# Code for all COX-regressions are the exact same just variables are changed.

# https://www.youtube.com/watch?v=arOdQrk09RI for general help with lifelines package!

import matplotlib
matplotlib.use('Agg')  # To save files
import matplotlib.pyplot as plt
import pandas as pd
from lifelines import KaplanMeierFitter, CoxPHFitter
import json

oligo_clean = pd.read_csv('data/TCGA/oligo_clean_survival.csv') # reload my prior csv fom data_loading_01
print(f"{len(oligo_clean)} patients") # making sure it is all good, 197 patients


# Onto major cox-regression for firstly clinical variables: https://numiqo.com/tutorial/cox-regression, https://lifelines.readthedocs.io/en/latest/Survival%20Regression.html for learning

cox_df = oligo_clean[['PFI.time', 'PFI', 
                 'age_at_initial_pathologic_diagnosis',
                 'neoplasm_histologic_grade',
                 'gender',
                 'radiation_therapy']].copy()

cox_df['PFI.time'] = cox_df['PFI.time'] / 365 # for years

# Label the variables as numbers (0,1)
cox_df['grade'] = (cox_df['neoplasm_histologic_grade'] == 'G3').astype(int)  # G3=1, G2=0
cox_df['male'] = (cox_df['gender'] == 'MALE').astype(int)  # Male=1, Female=0
cox_df['radiation'] = (cox_df['radiation_therapy'] == 'YES').astype(int)  # Yes=1, No=0

cox_df = cox_df.dropna(subset=['PFI.time', 'PFI', 'age_at_initial_pathologic_diagnosis']) #removing items with missing data points
print(f"Patient check: {len(cox_df)}") # hoping it stays at 197 patients (it does!!)

cox_input = cox_df[['PFI.time', 'PFI', 
                     'age_at_initial_pathologic_diagnosis',
                     'grade', 'male', 'radiation']].copy()


# Time to fit cox model (and understand the data slowly)
cph = CoxPHFitter()
cph.fit(cox_input, 
        duration_col='PFI.time', 
        event_col='PFI')

# Print results
cph.print_summary(decimals=3) # thats a lot of data... https://lifelines.readthedocs.io/en/latest/Survival%20Regression.html helped me understand it
# We see that gender does not really matter, radiation is unable to tell with the data, grade is significant, and age is certainly the most accurate prediction in terms of this batch of clinical data
# Great! Now we move on to the actual major biological signatures from my research:

# Load the mutations index: 
mutations = pd.read_csv('data/TCGA/LGG_mc3_mutations.txt.gz', 
                        sep='\t', 
                        index_col=0,
                        compression='gzip')

# Extract key oligodendroglioma markers
key_genes = ['IDH1', 'IDH2', 'TERT', 'CIC', 'FUBP1', 'ATRX', 'NOTCH1', 'PIK3CA']
# IDH1/2 are extremely prominent, as well as TERT. The others I classify as tier 2

found = [g for g in key_genes if g in mutations.index]
missing = [g for g in key_genes if g not in mutations.index] 
print(f"Found: {found}")
print(f"Missing: {missing}")

# We found all of them, so thats very good!


mut_subset = mutations.loc[found].T # Swap all the rows and columns transpose matrix (gemini snippit on how to transpose)
print(f"\nMutation rates in full data:") # which are most popular
print(mut_subset.mean().round(3))

# Merge using the clinical sampleID and the mutation index
oligo_combined = oligo_clean.merge(mut_subset, left_on='sampleID', right_index=True, how='inner')

print(f"Patients with both clinical and mutation data: {len(oligo_combined)}") # To merge them and see how many have it - 187



# Exact same COX approach for the mutations themselves:
# Prepare data for Cox regression
cox_data = oligo_combined[['PFI.time', 'PFI', 'IDH1', 'IDH2', 
                        'CIC', 'FUBP1', 'ATRX', 'NOTCH1', 
                        'PIK3CA']].copy()

# Drop any rows with missing values
cox_data = cox_data.dropna()
print(f"Patients for Cox regression: {len(cox_data)}") # Stays at the right amount

# Fit Cox model with IDH1 only first
cph = CoxPHFitter()
cph.fit(cox_data[['PFI.time', 'PFI', 'IDH1']], 
        duration_col='PFI.time', 
        event_col='PFI')

cph.print_summary()

# Tested it to work, so with IDH1 we see that those with it have a 57% reduction to risk of death, which aligns with research done in the field 
# The above was a mistake, it was increased for all GLIOMAS but for oligodendroglioma it shouldn't be 57%.

# Now, lets try it with all the genes I wrote down above

genes_to_test = ['IDH1', 'IDH2', 'TERT', 'CIC', 'FUBP1', 'ATRX', 'NOTCH1', 'PIK3CA'] # key_genes
univariable_results = []
print("Now running COX regression for all the major biological signatures!!!")

for gene in genes_to_test:
    test_df = oligo_combined[['PFI.time', 'PFI', gene]].dropna() # Incase one is missing the other genes, it doesnt throw out the entire dataset: https://lifelines.readthedocs.io/en/latest/Survival%20Regression.html helped me immensely for this first COX approach.
    test_df = test_df[test_df['PFI.time'] > 0] # As PFI can be =0
    test_df = test_df.copy()
    test_df['PFI.time'] = test_df['PFI.time'] / 365

    if test_df[gene].nunique() < 2:
        print(f"Skipping {gene}: No variation (either 100% mutated or 100% wildtype).") # Rule in COX regression and statistics. 
        continue

    try: # Time to run the actual cox fit on all the individual genes:
        cph_temp = CoxPHFitter()
        cph_temp.fit(test_df, duration_col='PFI.time', event_col='PFI')
        
        summary = cph_temp.summary
        hr = summary['exp(coef)'].values[0]         # Looking for both the HR (hazard risk) of the gene (which tells us what the gene exactly does to the survival risk)
        p_val = summary['p'].values[0]              # and the p value, the probability/how much we can 'trust' it. Lower = better
        
        univariable_results.append({
            'Gene': gene,
            'Hazard Ratio (HR)': round(hr, 3),
            'p-value': round(p_val, 4)                # Rounding these long decimal points for p and HR. Doing this individually to create a table
        })
    except Exception as e:
        print(f"Could not fit model for {gene}: {e}") # Incase of error 

results_df = pd.DataFrame(univariable_results).sort_values(by='p-value') # Sort the values by p-value to create a table letting me see the most prominent key-genes. 
print("\nGenomic Screening Resulst!!!")
print(results_df.to_string(index=False)) 

# We see that CIC is extremely prominent- 56% reduced risk of death
# IDH1 is also prominent as well
# IDH2 is also prominent!!
# for TERT, we see that the numbers are completely off- THIS IS BECAUSE:::: TERT mutations happen in the non-coding region, so it was not inside the dataset in the first place. This result cannot be trusted. 
# for TERT, we will have to use a seperate approach later.
# for the others, they are insignificant when looked at individually (either nothing individually, too minor, or small data set)



# Full model now incorporating both clinical features and IDH1 and CIC


# Map out our clinical features cleanly on the merged dataset
oligo_combined['grade'] = (oligo_combined['neoplasm_histologic_grade'] == 'G3').astype(int) 

# Now we set our variables as the clinical features that showed promise, as well as the IDH1 and CIC that also showed promise.
final_features = [
    'PFI.time', 'PFI', 
    'age_at_initial_pathologic_diagnosis', 'grade',
    'IDH1', 'CIC'
]

# Ensure we drop rows missing any of these specific features and convert time to years
final_cox_input = oligo_combined[final_features].dropna().copy()
cox_data = cox_data[cox_data['PFI.time'] > 0]
final_cox_input['PFI.time'] = final_cox_input['PFI.time'] / 365

print(f"Patients included in the final model: {len(final_cox_input)}")

# Fit it just like how we do it for all theo thers we did so far, final step!!!
cph_final = CoxPHFitter()
cph_final.fit(final_cox_input, duration_col='PFI.time', event_col='PFI')

print("\nFinal COX model results...")
cph_final.print_summary(decimals=3)

# We got 0.809 concordance which is super super good.
# It also says that CIC is more important than IDH1. Typically, CIC and IDH1 always happen together, so the fact that we were able to differentiate one from the other which such high concordance is superb.



# However, I also noticed a problem: the Oligodendroglioma might be misclassified, as relooking at previous numbers we see that IDH1 mutation was not there 100% of the time. I was mixed up in my knowledge. All Lower Grade Glioma has around that 'more-deathly' if it has IDH1, but I should be looking at Oligodendroglioma specifically, so it shouldn't even be a risk-factor. It should just be 100%.
# In this TCGA dataset, the scientists used a histological approach through a microscope, and basically they mislabelled them
# In order for something to be classified as Oligodendroglioma, it must have the IDH1/2 mutation.
# I decided to solve this using open-source software, cBioPortal, and clean my entire dataset and redo the cox regression.

# Redoing COX regression (essentially just copy and paste the prior cox with a different dataset)


# Time for cBioPortal quality control ;
with open('data/TCGA/lgggbm_sample_clinical.json', 'r') as f: # downloaded through cBioPortal website opensource
    sample_clinical = json.load(f) # https://www.geeksforgeeks.org/python/how-to-read-json-files-with-pandas/ for help with using pd. and json filesfrom cBioPortal
sample_df = pd.DataFrame(sample_clinical)

sample_wide = sample_df.pivot_table( # Pandas says we need to expand the matrix, so thats waht I did (help from gemini snippit to figure out why I was erroring)
    index='sampleId',
    columns='clinicalAttributeId',
    values='value',
    aggfunc='first'
).reset_index()  

mol_cols = ['sampleId', 'IDH_CODEL_SUBTYPE', 'IDH_STATUS', 'GRADE']
mol_data = sample_wide[mol_cols].copy() # We want to search ONLY for IDH mutation to limit it. Again, we need them ALL to have IDH to clean the set.


oligo_validated = oligo_combined.merge(mol_data,
                                       left_on='sampleID',
                                       right_on='sampleId',
                                       how='left') # Merging together the datasets. 

print(f"Raw cohort size: {len(oligo_validated)}")
print(oligo_validated['IDH_CODEL_SUBTYPE'].value_counts(dropna=False))
oligo_confirmed = oligo_validated[oligo_validated['IDH_CODEL_SUBTYPE'] == 'IDHmut-codel'].copy() #The filter itself

print(f"Confirmed Oligodendrogliomas: {len(oligo_confirmed)}")
print(f"Removed as gliomas: {len(oligo_validated) - len(oligo_confirmed)}")

# Now we should have a dataset of 127 Oligodendroglioma ONLY patients, the other 60 were misclassified (as before, they often misclassified them) which is honestly a pretty interesting finding (although accidental)!

# Onto basic COX, which is just copy-pasted from above with a correct dataset.
confirmed_univariable_results = []

for gene in genes_to_test:
    # Final clean and filtered dataset, oligo_confirmed
    test_df = oligo_confirmed[['PFI.time', 'PFI', gene]].dropna()
    test_df = test_df[test_df['PFI.time'] > 0] # As PFI = 0
    test_df = test_df.copy()
    test_df['PFI.time'] = test_df['PFI.time'] / 365


    if test_df[gene].nunique() < 2: # Exact same as before
        continue

    try:
        cph_temp = CoxPHFitter()
        cph_temp.fit(test_df, duration_col='PFI.time', event_col='PFI')
        
        summary = cph_temp.summary
        hr = summary['exp(coef)'].values[0]         
        p_val = summary['p'].values[0]              
        
        confirmed_univariable_results.append({
            'Gene': gene,
            'Hazard Ratio (HR)': round(hr, 3),
            'p-value': round(p_val, 4)                
        })
    except Exception as e:
        print(f"Could not fit model for {gene}: {e}")

confirmed_results_df = pd.DataFrame(confirmed_univariable_results).sort_values(by='p-value') 
print("\nGenomic Results (True Oligodendroglioma Cohort Only)!!!")
print(confirmed_results_df.to_string(index=False))


# Now that we have those out of the way, we are now going to see if these mutations are dependant on age and grade of the tumor like we did before. (but with a clean dataset...)
# We see that CIC is much further down the list, as it being tier 2 and actually not as prominent as IDH1! It just typically appears with IDH1, so that is unable to be used anymore.
# Also, we see that NOTCH1 is the HIGHEST on the list- this means that NOTCH1 actually is extremely prominent in this, if you have it you have double the risk of death!



# Same as what we did before
oligo_confirmed['grade_encoded'] = (oligo_confirmed['neoplasm_histologic_grade'] == 'G3').astype(int)

# And now we readd the important molecular signatures from the previous cox result!

final_strict_features = [
    'PFI.time', 'PFI', 
    'age_at_initial_pathologic_diagnosis', 'grade_encoded',
    'CIC', 'NOTCH1' # Added Notch1 and CIC as they were at the top of my list earlier!
]
strict_cox_input = oligo_confirmed[final_strict_features].dropna().copy()
cox_data = cox_data[cox_data['PFI.time'] > 0]
strict_cox_input['PFI.time'] = strict_cox_input['PFI.time'] / 365

cph_strict = CoxPHFitter()
cph_strict.fit(strict_cox_input, duration_col='PFI.time', event_col='PFI')

print("\nFinal model summary including age, grade, CIC, and NOTCH1")
cph_strict.print_summary(decimals=3)

# We see that age and grade, once again, dominate the landscape. Thus, second (tier 2) mutations matter but not as much as grade and age (of course)

# Really quickly saving dataset to do rna_analysis on the actual confirmed dataset 
oligo_confirmed.to_csv('data/TCGA/oligo_confirmed_master.csv', index=False)


# Now, lets graph the difference between grade 2 and 3 to see jsut how much it matters for survival:

fig, ax = plt.subplots(figsize=(10, 6)) #Same as data_loading_01.py setup! https://lifelines.readthedocs.io/en/latest/fitters/univariate/KaplanMeierFitter.html

# https://lifelines.readthedocs.io/en/latest/fitters/univariate/KaplanMeierFitter.html

# Filter our clean cohort into Grade 2 and Grade 3 buckets
grade2_mask = oligo_confirmed['neoplasm_histologic_grade'] == 'G2'  # We want to classify and seperate grade 2 and 3 from each other
grade3_mask = oligo_confirmed['neoplasm_histologic_grade'] == 'G3'

# Fit and plot Grade 2 Patients
kmf_g2 = KaplanMeierFitter()
kmf_g2.fit(
    durations=oligo_confirmed[grade2_mask]['PFI.time'] / 365,  # Convert days to years
    event_observed=oligo_confirmed[grade2_mask]['PFI'],
    label=f"WHO Grade II (n={grade2_mask.sum()})"
)
kmf_g2.plot_survival_function(ax=ax, color='#2c7bb6', ci_show=True)

# Fit and plot Grade 3 Patients
kmf_g3 = KaplanMeierFitter()
kmf_g3.fit(
    durations=oligo_confirmed[grade3_mask]['PFI.time'] / 365, # days to years again 
    event_observed=oligo_confirmed[grade3_mask]['PFI'],
    label=f"WHO Grade III (n={grade3_mask.sum()})"
)
kmf_g3.plot_survival_function(ax=ax, color='#d7191c', ci_show=True)

# 3. Basic formats and image boundaries
ax.set_xlabel('Time (Years)', fontsize=12, fontweight='bold')
ax.set_ylabel('Progression-Free Survival Probability', fontsize=12, fontweight='bold')
ax.set_title('Progression-Free Interval (PFI) by Histologic Grade\nin Molecularly Confirmed Oligodendroglioma ($IDH$mut-codel)',
             fontsize=14, fontweight='bold', pad=15)

# 50% survival baseline (for median and image details)
ax.axhline(y=0.5, color='gray', linestyle='--', alpha=0.5)

# axes limits + grid
ax.set_ylim(0, 1.05)
ax.set_xlim(0, oligo_confirmed['PFI.time'].max() / 365)
ax.grid(True, linestyle=':', alpha=0.6)
ax.legend(fontsize=11, loc='lower left')

plt.tight_layout()

# Saves
plt.savefig('figures/02_TCGA_confirmed_grade_KM.png', dpi=300)


# Currently at home right now! In the rather long bus ride home which took upwards of an hour, I realized that...
# I was using OS time (which is essentially time of death), instead of PFI: The industry standard 'Progression-Free Interval'
# This tracks essentially the length they have had the tumor, not the age of the patient itself. The reason why age was always p<0.005 was because one of the variables was both the indepedent and the dependant!
# With PFI, I see that age still is very prominent, but NOTCH1 is also SUPER SUPER important as well, a really novel finding!
# This stays in the independent variable as well as the mutivariable cox models. Now I will validate this in CCGA findings on my next python file.


# https://www.youtube.com/watch?v=arOdQrk09RI at timestamp 11:09 taught me about Shoenfeld residuals, and gave me the check if NOTCH1 or anything is time-dependant!
# https://www.youtube.com/watch?v=FQbNv0sZBYc https://www.google.com/url?sa=t&source=web&rct=j&opi=89978449&url=https://stats.stackexchange.com/questions/547078/schoenfeld-residuals-plain-english-explanation-please&ved=2ahUKEwiDudTey8GVAxWHBTQIHQ0qI_QQFnoECCIQAQ&usg=AOvVaw2daKw8m8VUePusNgiPA4q9
# Essentially, the schoenfeld residues are basically testing if the hazard risks / risk effects stay the entire time. Making sure they aren't time dependant in a sense.
# I do not really understand the specifics too well, but that is the general idea.
# https://timeseriesreasoning.com/contents/schoenfeld-residuals/ another great source

from lifelines.statistics import proportional_hazard_test #https://lifelines.readthedocs.io/en/latest/jupyter_notebooks/Proportional%20hazard%20assumption.html helped a lot
# High p-value >0.05 = assumption stays (good)
# Low p-value <0.05 = effect changes over time, need time-varying model (results inaccurate)
results_ph = proportional_hazard_test(cph_strict,
                                      strict_cox_input,
                                      time_transform='rank')

print("\nPH Test (Schoenfeld Residuals)")
print(results_ph.summary)
# P values of 0.79, 0.52, 0.08, and 0.95 (all >0.05)! This means that we are a-okay for these results. Just an additional expansion I was interested in and stumbled upon.




