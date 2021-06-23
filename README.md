# bulk-rnaseq
Simple  workflow to quantify gene-level RNA abundance and detect differentially expressed genes (DEGs) 
from bulk RNAseq samples. The pipeline uses `kallisto` or `salmon` to quantify transcript level abundance and `DESeq2` 
to normalize counts and detect DEGs. 

## Installation
1. Make sure `conda`, `mamba`, and `snakemake` are installed
2. Download the appropriate [`kallisto`](https://github.com/pachterlab/kallisto-transcriptome-indices/releases) or  [`salmon`](http://refgenomes.databio.org) references or build your own
3. Clone this repository or a create new repository from this template
4. Describe your samples in `samples.csv`
5. Modify the settings in `config.yaml`
6. (Optional) If you plan on using a SLURM cluster, fill out the `#SBATCH` directives in `run_pipeline.sh` and the `out` and `account` fields in `cluster.json`. The `out` path must already exist; Snakemake will not create directories for you
8. (Optional) If you want to run the pipeline in a Singularity environment for full reproducibility, install Singularity. 
`run_pipeline.sh` assumes Singularity is installed. Delete the `--use-singularity` flag if you want to skip using a Singularity environment.
The same goes for using `conda` environments, although it is recommended to use both for the best reproducibility.

### conda vs mamba
`snakemake` recommends the `mamba` installer instead of `conda` because of `mamba`'s better performance. See the [documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) for installing `mamba`. If you prefer to use `conda`, see the `snakemake`
command line arguments to install environments using `conda`.

## Salmon vs Kallisto
I recommend using `salmon` for quantification as it automatically infers the library strandedness and supports automatic gene annotation via `tximeta`.
If you prefer to use `kallisto`, make sure to specify the strandedness in the sample file.

## Sample File Configuration
`bulk-rnaseq` requires you specify a CSV file with a metadata about samples to analyze. 
There are several required columns:
* `patient`
* `condition`
* `fq1`
* `fq2`
* `strandedness` (required by `kallisto`; optional if using `salmon`)

The `patient` and `condition` columns are used to create a sample specific ID for downstream processing and identification.
The values specified in the these columns should create unique values when combined as `{patient}-{condition}`. 
`condition` should specify the experimental condition of the sample (`normal` or `tumor`, `primary` or `metastatic` etc). 

**NOTE** `salmon` will automatically infer the library stranding and choose the best option. `kallisto` requires you specify
the stranding info. See [`check-strand`](https://github.com/tjbencomo/check-strand) to determine the proper stranding if you prefer to
run `kallisto`.


## Running the pipeline
Run jobs locally
```
snakemake -j [cores] --use-singularity --use-conda
```
Run jobs on a SLURM cluster
```
sbatch run_pipeline.sh [path to directory with Snakefile]
```

**NOTE:** It is highly recommended you run `snakemake` with the `--use-conda` and `--use-singularity` flags for
maximum reproducibility

## Output
Results are stored in `results/` and include:
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
