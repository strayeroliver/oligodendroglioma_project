# Oligodendroglioma Survival Analysis
## Description
The purpose of this project is to analyze the differences in oligodendroglioma patients to understand why some live for 4 months and others can survive up to 20 years.  

The end goal is to analyze:  
1. Clinical details about the patient (Age, grade of tumor, chemotherapy/radiation)
2. Major biological signatures of the patient (Well-known like IDH1/2, TERT, etc.)
3. Full RNA-sequence of the patient (Entire genome sequence to find new, novel data never done before)
4. Ratio of tumor suppressors vs Secondary Oncogenes 
5. The pathways behind the genes in steps 2/3 to explain our answers
6. Methylation data of all patients, describing their position near the promoter region

To find new, novel answers and create a multivariable tool that clinicians can use to predict oligodendroglioma patient's survival length and prescribe the correct treatment.
 
## Prerequisites/Setup

Create and activate a virtual environment:  
  ```
   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
   bash Miniconda3-latest-Linux-x86_64.sh
```
Then restart your teminal:
```
   conda create -n oligodendroglioma python=3.11
   conda activate oligodendroglioma
   ```

Install required python packages:  
```
pip install numpy pandas scipy matplotlib seaborn scikit-learn lifelines gseapy
```

You will also need R for this project, so:
```
sudo apt update
sudo apt install r-base r-base-dev -y
```

## Clone the Repository
To get the code and data required:
```
git clone https://github.com/strayeroliver/oligodendroglioma_project.git
cd oligodendroglioma_project
```
## Running the Analysis
To follow along with my journey, simply run the notebooks in the following order:

```
python notebooks/data_loading_01.py
python notebooks/cox_regression_02.py
python notebooks/CGGA_validatoin_03.py
python notebooks/rna_analysis_04.py
Rscript notebooks/LASSO_cox_04.R 
python notebooks/...
```

And so forth.  
Rscripts are only needed if you would like to see the process of doing LASSO-cox and how I got my values.

