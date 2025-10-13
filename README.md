This is the repository for the metadata, sequencing/proteomic data, and code necessary for data analysis and figure generation for the manuscript titled "Distinct Host and Microbial Biology Distinguishes Lower Respiratory Tract Infection from Incidental Pathogen Carriage". Link to manuscript will be added once accepted/published.

#File descriptions
1. metadata.csv: includes limited, de-identified metadata for the entire cohort, including age, sex, clinical adjudication, clinical microbiology, mNGS hits (and the sum NT-RPM of the pathogenic viruses and bacteria detected), and shannon diversity index.
2. proteomics_data.csv: Relative FABP4 protein value for each participant (if available), as measured using the SomaScan 7k assay.
3. gene_counts.csv: Gene counts table, filtered to include genes that had at least 10 counts in 20% of the samples in the cohort. 

#Code descriptions
1. host_analyses_code.R: R script for performing all of the host- and integrated host-microbial analyses.
2. classifier_code.R: R script for generating host, microbial, and integrated classifiers.
