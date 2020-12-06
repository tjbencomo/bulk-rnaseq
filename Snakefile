singularity: "docker://continuumio/miniconda3"

import sys
import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version

configfile: "config.yaml"
samples_fp = config['samples']
samples = pd.read_csv(samples_fp, dtype=str)
samples['id'] = samples['patient'] + '-' + samples['sample']
samples = samples.set_index(["id"], drop=False)
samples = samples.sort_index()

kallisto_idx = config['kallisto_index']
kallisto_tx2g = config['kallisto_tx2gene']
kallisto_threads = config['kallisto_threads']

# parse factor levels into readable string format
# pass this string as param to deseq2_init.R
levels = config['factor_levels']
lvlstr = []
for level in levels:
    lvls = levels[level].replace(" ", "")
    lvlstr.append(f"{level}={lvls}")
var_levels = ';'.join(lvlstr)

def get_fqs(wildcards):
    if isPE(wildcards):
        return {
            'fq1' : samples.loc[(wildcards.sample_id), 'fq1'],
            'fq2' : samples.loc[(wildcards.sample_id), 'fq2']
            }
    else:
        return {'fq1' : samples.loc[(wildcards.sample_id), 'fq1']}

def getStrand(wildcards):
    s = samples.loc[wildcards.sample_id, 'strandedness'].lower()
    if s == 'forward':
        return '--fr-stranded'
    elif s == 'reverse':
        return '--rf-stranded'
    elif s == 'none' or s == 'unstranded':
        return ''
    else:
        raise ValueError(f"Unrecognized strand type for {wildcards.sample_id}")

def isPE(wildcards):
    if pd.isna(samples.loc[wildcards.sample_id, 'fq2']):
        return False
    else:
        return True

###########################################################
### Snakemake Rules
###########################################################

rule targets:
    input:
        expand("kallisto/{sample_id}", sample_id=samples['id']),
        #expand("results/{contrast}/{estimate}_foldchanges.csv", contrast=config['contrasts'], estimate=['mle', 'map']),
        #expand("results/{contrast}/ma_plot.svg", contrast=config['contrasts']),
        #"results/pca_plot.svg",
        #"results/normalized_counts.rds",
        "deseq2/all.rds"

rule kallisto:
    input:
        unpack(get_fqs),
        idx=kallisto_idx
    output:
        directory("kallisto/{sample_id}")
    params:
        strand = lambda wildcards: getStrand(wildcards),
        fqs = lambda wildcards, input: input.fq if not isPE(wildcards) else f"{input.fq1} {input.fq2}"
    log:
        "logs/kallisto/{sample_id}.log"
    threads: 4
    conda:
        "envs/quant.yml"
    shell:
        """
        kallisto quant -i {input.idx} -o {output} -t {threads} {params.strand} \
            {params.fqs} | tee {log}
        """

rule deseq2_init:
    input:
        cts = expand("kallisto/{sample_id}", sample_id = samples['id']),
        samples=samples_fp,
        tx2g = kallisto_tx2g
    output:
        deseq="deseq2/all.rds",
        cts="results/normalized_counts.rds"
    params:
        levels=var_levels
    script:
        "scripts/deseq2.R"
