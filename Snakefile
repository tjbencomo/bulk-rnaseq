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
se_frag_length = config['single_end_frag_length']
se_frag_sd = config['single_end_frag_sd']

# parse factor levels into readable string format
# pass this string as param to deseq2_init.R
levels = config['factor_levels']
lvlstr = []
for level in levels:
    lvls = levels[level].replace(" ", "")
    lvlstr.append(f"{level}={lvls}")
var_levels = ';'.join(lvlstr)

design_formula = config['design_formula']

pca_labels = config['pca_labels'].replace(" ", "")

def get_fqs(wildcards):
    if isPE(wildcards):
        return {
            'fq1' : samples.loc[(wildcards.sample_id), 'fq1'],
            'fq2' : samples.loc[(wildcards.sample_id), 'fq2']
            }
    else:
        return {'fq' : samples.loc[(wildcards.sample_id), 'fq1']}

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

def get_contrast(wildcards):
    return config['contrasts'][wildcards.contrast]

###########################################################
### Snakemake Rules
###########################################################

rule targets:
    input:
        expand("kallisto/{sample_id}", sample_id=samples['id']),
        expand("results/{contrast}/{estimate}_foldchanges.csv", contrast=config['contrasts'], estimate=['mle', 'map']),
        expand("results/{contrast}/{estimate}_ma.svg", contrast=config['contrasts'], estimate=['mle', 'map']),
        "results/pca_plot.svg",
        "results/normalized_counts.rds",
        "deseq2/all.rds",
        "qc/multiqc_report.html"

rule kallisto:
    input:
        unpack(get_fqs),
        idx=kallisto_idx
    output:
        directory("kallisto/{sample_id}")
    params:
        strand = lambda wildcards: getStrand(wildcards),
        fqs = lambda wildcards, input: input.fq if not isPE(wildcards) else f"{input.fq1} {input.fq2}",
        single = lambda wildcards: '' if isPE(wildcards) else '--single',
        frag_length = lambda wildcards: '' if isPE(wildcards) else f"-l {se_frag_length}",
        frag_sd = lambda wildcards: '' if isPE(wildcards) else f"-s {se_frag_sd}"
    log:
        "logs/kallisto/{sample_id}.log"
    threads: 4
    conda:
        "envs/quant.yml"
    shell:
        """
        kallisto quant -i {input.idx} -o {output} -t {threads} \
            {params.strand} {params.single} {params.frag_length} {params.frag_sd} \
            {params.fqs} &> >(tee {log})
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
        formula=design_formula,
        levels=var_levels
    log:
        "logs/deseq2/init.log"
    conda:
        "envs/deseq2.yml"
    script:
        "scripts/deseq2.R"

rule diffexp:
    input:
        "deseq2/all.rds"
    output:
        mleres = "results/{contrast}/mle_foldchanges.csv",
        mapres = "results/{contrast}/map_foldchanges.csv",
        mlema = "results/{contrast}/mle_ma.svg",
        mapma = "results/{contrast}/map_ma.svg"
    params:
        formula=design_formula,
        contrast = get_contrast
    log:
        "logs/deseq2/{contrast}.log"
    conda:
        "envs/deseq2.yml"
    script:
        "scripts/diffexp.R"

rule multiqc:
    input:
        expand("logs/kallisto/{sample_id}.log", sample_id=samples['id'])
    output:
        "qc/multiqc_report.html"
    log:
        "logs/multiqc.log"
    conda:
        "envs/qc.yml"
    wrapper:
        "0.50.4/bio/multiqc"

rule pca:
    input:
        "results/normalized_counts.rds"
    output:
        "results/pca_plot.svg"
    params:
        label_vars = pca_labels
    log:
        "logs/plot_pca.log"
    conda:
        "envs/deseq2.yml"
    script:
        "scripts/plot_pca.R"
