## Phenotype data
The phenotype data can be found on the galaxy server.

### Data dictionary
ftp://ftp.ncbi.nlm.nih.gov/dbgap/studies/phs000424/phs000424.v7.p2/pheno_variable_summaries/phs000424.v7.pht002742.v7.GTEx_Subject_Phenotypes.data_dict.xml

## R Scripts
R code can be found in `src` folder.

#### Deploying code to the Galaxy server
1. Commit the code and push to GitHub.
2. Run the following command on your local machine from the project directory:
```
scp -P 2222 -r src/ designproject@galaxy.ugent.be:/data/designproject/blood_clonal_expansion/
```
Do not run analysis if your code hasn't been pushed to GitHub, since this compromises reproducibility.

#### Running processing with MAGE
Run the command following from the `blood_clonal_expansion` folder on the server:
```
Rscript mage_processing.R mage_processing.yaml
```

`mage_processing.yaml` should have following structure:
```
counts_path: "/data/designproject/results/sequences/"
seqem_path: "/data/designproject/sequem/"
results_path: "/data/designproject/mage_results/"
chromosomes: [1, 2, 3, 4, X]
parallelism: 4
```

#### Running statistical analysis

Run the command following from the `blood_clonal_expansion` folder on the server:
```
Rscript statistical_analysis.R statistical_analysis.yaml
```

`statistical_analysis.yaml` should have following structure:
```
counts_paths: "/data/designproject/results_combined/sequences/"
phenotypes_path: "/data/designproject/Phenotypes.txt"
mage_results_path: "/data/designproject/mage_results/"
plots_path: "/data/designproject/stat_analysis_results/"
results_path: "/data/designproject/stat_analysis_results/"
chromosomes: ["1"]
paralellism: 4
```


