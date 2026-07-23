# Data Sources

This project does not commit TCGA or CGGA data to the repo (due to privacy laws for patients). Follow the steps below to download everything into data/TCGA/ and data/CGGA/.

## 1. TCGA data (from UCSC Xena)

Go to: https://xenabrowser.net/datapages/?cohort=TCGA%20Lower%20Grade%20Glioma%20(LGG)

Download these files and put them in data/TCGA/:

- Under Phenotype, download "survival/LGG_survival.txt". Save as LGG_survival.txt (it downloads as .gz, unzip this manually)
- Under Phenotype, download "LGG_clinicalMatrix". Save as LGG_clinicalMatrix (also unzip after downloading)
- Under gene expression RNAseq, download "HiSeqV2". Save as LGG_RNAseq.gz (keep it zipped)
- Under DNA methylation, download "HumanMethylation450". Save as LGG_methylation450k.gz (keep it zipped)
- Under somatic mutation (pan-cancer), search for "mc3" and download "LGG_mc3.txt". Save as LGG_mc3_mutations.txt.gz (keep it zipped)

## 2. Molecular subtype data (from cBioPortal)

Go to: https://www.cbioportal.org/study/summary?id=lgggbm_tcga_pub

Click the "Clinical Data" tab and export the sample-level clinical data as JSON. (its the little download button next to # of patients on the right side of the top bar)

Save as data/TCGA/lgggbm_sample_clinical.json

## 3. Illumina 450k manifest

Go to: https://support.illumina.com/downloads/infinium_humanmethylation450_product_files.html

Download "HumanMethylation450_15017482_v1-2.csv" (the manifest CSV, not the .bpm file).

Save as data/TCGA/humanmethylation450_15017482_v1-2.csv

## 4. CGGA data

Go to: http://www.cgga.org.cn/download.jsp

Under mRNAseq_325, download the "Clinical Data" and "Expression Data (RSEM)" files (dated 20200506):
- Save clinical file as data/CGGA/CGGA.mRNAseq_325_clinical.20200506.txt
- Save expression file as data/CGGA/CGGA.mRNAseq_325.RSEM-genes.20200506.txt

Under mRNAseq_693, do the same:
- Save clinical file as data/CGGA/CGGA.mRNAseq_693_clinical.20200506.txt
- Save expression file as data/CGGA/CGGA.mRNAseq_693.RSEM-genes.20200506.txt

Under WEseq_286, download the clinical file and the mutation file:
- Save clinical file as data/CGGA/CGGA.WEseq_286_clinical.20200506.txt
- Save mutation file as data/CGGA/CGGA.WEseq_286.20200506.txt

If any of them are zipped, please unzip these.
