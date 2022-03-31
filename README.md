# Introduction

TCRP (Transfer of Cellular Response Prediction) is a few-shot machine learning algorithm to perform transfer learning for drug response predictions in cancer. It is used to train a neural network model using data in one experimental context (e.g. cell lines) that can be readily tuned to new contexts (e.g. patients) using few additional samples. In pilot experiments, the model has been shown to quickly adapt when switching among different tissue types and in moving from cell-line models to clinical contexts, including patient-derived tumor cells and patient-derived xenografts.

The associated publication for this method can be found via this citation:

Ma J, Fong SH, Luo Y., Bakkenist CJ, Shen JP, Mourragui S, Wessels LFA, Hafner M, Sharan R, Peng J, Ideker T.  Few-shot learning creates predictive models of drug response that translate from high-throughput screens to individual patients. Nat Cancer. 2021 Feb;2(2):233-244. doi: 10.1038/s43018-020-00169-2. Epub 2021 Jan 25. PMID: 34223192 [Pubmed](https://pubmed.ncbi.nlm.nih.gov/34223192/)

This GitHub repository provides an implementation of data preprocessing to generate model input in an attempt to reproduce results from Challenge 2 of the publication: transfer to PDTCs.

# Running the code

Install the specified dependencies in requirements.txt, or run the Dockerfile in the environment folder. 

Generating input for the TCRP model is as simple as running the notebook "code/Process_merged.ipynb". This notebook will receive both GDSC and PDTC dataset files (provided) as input, and will generate drug-specific feature and label files under the path "/data/merged". To run the TCRP model with this data, you must take this resutling folder and move it to this [TCRP model repository](https://github.com/emilyso-99/TCRP_pipeline)(instructions on paths are outlined)


## Data Availability: 

The gene expression and somatic mutation profiles used in the paper for each cell line are from the Cancer Cell Line Encyclopedia (CCLE) project. These data can be downloaded from the DepMap website: http://depmap.org/portal/download/. 
The drug response data used in the paper for each cell line can be downloaded from the GDSC 1000 website: http://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources/.

Expression data, somatic mutations and drug responses for the analysis of Patient-Derived Tumor Cells (PDTC) reported in the paper can be downloaded from the following URL: http://figshare.com/articles/Bruna_et_al_A_biobank_of_breast_cancer_explants_with_preserved_intra-tumor_heterogeneity_to_screen_anticancer_compounds_Cell_2016/2069274. 

All relevant data for the Patient Derived Xenograft (PDX) models can be extracted from Supplementary Table 1 of the paper ‘High-throughput screening using patient-derived tumor xenografts to predict clinical trial drug response.’
