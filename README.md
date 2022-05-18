# bulk-rnaseq
Simple  workflow to quantify gene-level RNA abundance and detect differentially expressed genes (DEGs) 
from bulk RNAseq samples. The pipeline uses `kallisto` or `salmon` to quantify transcript level abundance and `DESeq2` 
to normalize counts and detect DEGs. 

## Resource Requirements
Salmon usually works with 8GB to 12GB of RAM. 16GB is recommended for DESeq2 with small datasets. Larger amounts of memory
will be required for larger datasets. 

## Installation
1. Install Singularity and `snakemake`
2. Download the [`salmon`](http://refgenomes.databio.org) references or build your own
3. Clone this repository or a create new repository from this template
4. Describe your samples in `samples.csv`
5. Modify the settings in `config.yaml`
6. (Optional) If you plan on using a SLURM cluster, fill out the `#SBATCH` directives in `run_pipeline.sh` and the `account` field in `cluster.json`. 

## Salmon vs Kallisto
As of April 2022, `bulk-rnaseq` only supports Salmon quantification. Older versions supported Kallisto,
but this was dropped due to better Salmon integration with tximeta

## Sample File Configuration
`bulk-rnaseq` requires you specify a CSV file with a metadata about samples to analyze. 
There are several required columns:
* `patient`
* `condition`
* `fq1`
* `fq2`

The `patient` and `condition` columns are used to create a sample specific ID for downstream processing and identification.
The values specified in the these columns should create unique values when combined as `{patient}-{condition}`. 
`condition` should specify the experimental condition of the sample (`normal` or `tumor`, `primary` or `metastatic` etc). 

**NOTE** `salmon` will automatically infer the library stranding and choose the best option.

## Running the pipeline
Run jobs locally
```
snakemake -j [cores] --use-singularity
```
Run jobs on a SLURM cluster
```
# Must run within the bulk-rnaseq directory
sbatch run_pipeline.sh
```

## Output
Results are stored in `results/[gene|transcript]-level` and include:
* `normalized_counts.rds` - DESeq2 object with VST normalized log2 counts. Useful for downstream visualization, clustering, and machine learning
* `pca_plot.svg` - Plot of the samples described by the first two principal components. Useful to check for batch effects and identify bad samples

Each contrast has its own folder in `results/` with the following files:
* `mle_foldchanges.csv` - List of log2 fold changes estimated by DESeq2 using maximum likelihood estimation (MLE)
* `map_foldchanges.csv` - List of log2 fold changes estimated by DESeq2 via shrinkage using the `apeglm` package
* `mle_ma.svg` - MA plot showing the MLE fold changes
* `map_ma.svg` - MA plot showing the MAP fold changes

The full DESeq2 object used for differential expression analysis can be found in `deseq2/`.

A MultiQC html report with quantification statistics and FASTQC reports for each sample is located in `qc/`. 
Check this to verify a reasonable proportion of reads were pseudoaligned.
`validateFastq_summary.csv` has the results of running `qc/validateFastq` on each sample's input FASTQ files.
This checks for proper FASTQ formatting and properly sorted read headers.
