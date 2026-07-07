# For data importation and basic understanding of said data


import matplotlib
matplotlib.use('Agg')  # Non-interactive backend - saves files without displaying
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter # Had to do research on how to use Kaplan Meier Fitter https://lifelines.readthedocs.io/en/latest/fitters/univariate/KaplanMeierFitter.html

import pandas as pd
import os
os.chdir('/home/oliver/oligodendroglioma_project')

#loading in and analysis of TCGA dataset: https://pmc.ncbi.nlm.nih.gov/articles/PMC6066282/
# https://xenabrowser.net/datapages/?cohort=TCGA%20Lower%20Grade%20Glioma%20(LGG)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443

survival = pd.read_csv('data/TCGA/LGG_survival.txt', sep='\t') #load data in

# Actual qualities:

print("dataset qualities survival")
print(f"Total Patients: {len(survival)}") # output 529 total LGG patients, I want Oligodendroglioma specific
print(f"Patients Died (OS = 1): {survival['OS'].sum()}") # number of deceased patients- good for dataset, OS = 1 if dead, 0 if alive
print(f"Patients Alive: {(survival['OS']==0).sum()}") #529-133 = 396, just checking alive
print(f"Survival time (days):")
print(survival['OS.time'].describe()) #avg of 984 days, max is 6423 days (20 years)
print(f"Survival time (years)")
print((survival['OS.time']/365).describe()) #2.7 avg, 17.6 max to be precise

# Now we match up each patient with clinical data:

clinical = pd.read_csv('data/TCGA/LGG_clinicalMatrix', sep='\t') #load it in
print(f"Shape: {clinical.shape}")
print("\nHistology Types:")
print(clinical['histological_type'].value_counts()) # total 198 oligodendroglioma
oligo = clinical[clinical['histological_type'] == 'Oligodendroglioma'].copy()
print(f"Number of Oligodendroglioma Patients: {len(oligo)}") #198

#Now merge the clinical data onto this survival length
oligo_survival = oligo.merge(survival, left_on='sampleID', right_on='sample', how='inner') #Merge came from: https://jakevdp.github.io/PythonDataScienceHandbook/03.07-merge-and-join.html guide. Understood and carried it out.
print(f"After merging with survival data: {len(oligo_survival)}")  #Checking to make sure, 198 patients

print("\nOLIGODENDROGLIOMA SURVIVAL (years):")
print((oligo_survival['OS.time']/365).describe()) #avg 3 years, max 15.19, oligodendroglioma only

print(f"\nDied: {int(oligo_survival['OS'].sum())}")
print(f"Alive: {int((oligo_survival['OS']==0).sum())}") #47 deaths, 151 alive

# remove the patients with underfined patient lifespan (essentially alive)
oligo_clean = oligo_survival.dropna(subset=['OS.time', 'OS']).copy() #as NaN for time is inside patient data
print(f"Patients after removing missing survival data: {len(oligo_clean)}")




# Time for visualization:
# https://projector-video-pdf-converter.datacamp.com/25923/chapter2.pdf helped vastly
kmf = KaplanMeierFitter() # https://lifelines.readthedocs.io/en/latest/fitters/univariate/KaplanMeierFitter.html
kmf.fit(
    durations=oligo_clean['OS.time']/365, #in years
    event_observed=oligo_clean['OS'],  #Survival time
    label='Oligodendroglioma (n=197)'
)

fig, ax = plt.subplots(figsize=(10, 6)) # size of graph
kmf.plot_survival_function(ax=ax, ci_show=True) # quick gemini snippet (1 line) to plotting it


#Set my graph labels
# Heavily relied on https://lifelines.readthedocs.io/en/latest/fitters/univariate/KaplanMeierFitter.html for guidance in this, as well as https://projector-video-pdf-converter.datacamp.com/25923/chapter2.pdf

ax.set_xlabel('Time (Years)', fontsize=12)
ax.set_ylabel('Survival Probability', fontsize=12)
ax.set_title('Overall Survival -- TCGA Oligodendroglioma', fontsize=14)
ax.axhline(y=0.5, color='red', linestyle='--', alpha=0.5, label='50% survival') # for median survival time, line at 50%
ax.legend(fontsize=11)
plt.tight_layout()
plt.savefig('figures/01_overall_survival_KM.png', dpi=150) #saving it -- graph turned out very good
plt.show()

print(f"Median survival: {kmf.median_survival_time_:.2f} years") #finding median survival, 7.96 years https://projector-video-pdf-converter.datacamp.com/25923/chapter2.pdf

print(oligo_clean['neoplasm_histologic_grade'].value_counts(dropna=False)) # Looking at actual features that could contribute to survival length prior to looking at individual biological signatures/rna methyylation analysis
oligo_clean.to_csv('data/TCGA/oligo_clean_survival.csv', index=False) 
print(f"Saved {len(oligo_clean)} patients to data/TCGA/oligo_clean_survival.csv") # https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.value_counts.html









