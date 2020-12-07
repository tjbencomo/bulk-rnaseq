# bulk-rnaseq
Simple  workflow to quantify gene-level RNA abundance and detect differentially expressed genes (DEGs) 
from bulk RNAseq samples. The pipeline uses `kallisto` to quantify transcript level abundance and `DESeq2` 
to normalize counts and detect DEGs. 

## Installation
1. Install Anaconda or Miniconda and then `conda install snakemake`
2. Download the appropriate kallisto references from [here](https://github.com/pachterlab/kallisto-transcriptome-indices/releases) or build your own
3. Clone the repository
4. Describe your samples in `samples.csv`
5. Modify the settings in `config.yaml`
6. (Optional) If you plan on using a SLURM cluster, fill out the `#SBATCH` directives in `run_pipeline.sh` and the `out` field in `cluster.json`
7. (Optional) If you want to run the pipeline in a Singularity environment for full reproducibility, install Singularity. 
`run_pipeline` assumes Singularity is installed. Delete the `--use-singularity` flag if you want to skip using a Singularity environment.

## Running the pipeline
Run jobs locally
```
snakemake -j [cores] --use-singularity --use-conda
```
Run jobs on a SLURM cluster
```
sbatch run_pipeline.sh
```
