singularity: "docker://condaforge/mambaforge"
# singularity: "docker://continuumio/miniconda3"

import sys
import pandas as pd
from pathlib import Path
from snakemake.utils import validate
from snakemake.utils import min_version

configfile: "config.yaml"
samples_fp = config['samples']
samples = pd.read_csv(samples_fp, dtype=str)
samples['id'] = samples['patient'] + '-' + samples['condition']
samples = samples.set_index(["id"], drop=False)
samples = samples.sort_index()

# Logs
slurm_logdir = config['slurm_log_dir']
logpath = Path(slurm_logdir)
logpath.mkdir(parents=True, exist_ok=True) 

# quant_program = config['quant_program']
quant_program = 'salmon'

# kallisto_idx = config['kallisto_index']
qthreads = config['quant_threads']
dthreads = config['deseq_threads']
# se_frag_length = config['single_end_frag_length']
# se_frag_sd = config['single_end_frag_sd']
salmon_idx = config['salmon_index']
salmon_add_reps = config['salmon_compute_replicates']

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

salmon_env = config['salmon_container']
# kallisto_env = config['kallisto_container']
fastqc_env = config['fastqc_container']
multiqc_env = config['multiqc_container']
r_env = config['r_container']
validateFastq_env = config['validateFastq_container']

def get_fqs(wildcards):
    if isPE(wildcards):
        return {
            'fq1' : samples.loc[(wildcards.sample_id), 'fq1'],
            'fq2' : samples.loc[(wildcards.sample_id), 'fq2']
            }
    else:
        return {'fq' : samples.loc[(wildcards.sample_id), 'fq1']}

# def getStrand(wildcards):
#     if quant_program == 'salmon':
#         return 'Not using kallisto'
#     else:
#         if 'strandedness' not in samples.columns:
#             raise ValueError("Set to use kallisto for quantification but not stranding specified!")
#     s = samples.loc[wildcards.sample_id, 'strandedness'].lower()
#     if s == 'forward' or s == 'stranded':
#         return '--fr-stranded'
#     elif s == 'reverse':
#         return '--rf-stranded'
#     elif s == 'none' or s == 'unstranded':
#         return ''
#     else:
#         raise ValueError(f"Unrecognized strand type for {wildcards.sample_id}")

def isPE(wildcards):
    if pd.isna(samples.loc[wildcards.sample_id, 'fq2']):
        return False
    else:
        return True

def get_quants(wildcards):
    # if quant_program == 'kallisto':
        # quants = [f"kallisto/{sid}" for sid in samples['id']]
    if quant_program == 'salmon':
        quants = [f"salmon/{sid}" for sid in samples['id']]
    else:
        raise ValueError("quant_program must be 'salmon'")
    return {'cts' : quants}

def get_contrast(wildcards):
    return config['contrasts'][wildcards.contrast]

def get_multiqc_input(wildcards):
    # if quant_program == 'kallisto':
        # logs = [f"logs/kallisto/{sample_id}.log" for sample_id in samples['id']]
    # else:
    logs = [f"salmon/{sample_id}" for sample_id in samples['id']]
    fastqc = [f"qc/fastqc/{sample_id}" for sample_id in samples['id']]
    return logs + fastqc


###########################################################
### Snakemake Rules
###########################################################

rule targets:
    input:
        expand("{program}/{sample_id}", program=quant_program, sample_id=samples['id']),
        expand("results/{contrast}/{estimate}_foldchanges.csv", contrast=config['contrasts'], estimate=['mle', 'map']),
        expand("results/{contrast}/{estimate}_ma.svg", contrast=config['contrasts'], estimate=['mle', 'map']),
        expand("qc/validateFastq/{sample_id}.txt", sample_id=samples['id']),
        "results/pca_plot.svg",
        "results/normalized_counts.rds",
        "deseq2/all.rds",
        "deseq2/se_tx.rds",
        "qc/multiqc_report.html",
        "qc/validateFastq_summary.csv"

# rule kallisto:
#     input:
#         unpack(get_fqs),
#         idx=kallisto_idx
#     output:
#         directory("kallisto/{sample_id}")
#     params:
#         strand = lambda wildcards: getStrand(wildcards),
#         fqs = lambda wildcards, input: input.fq if not isPE(wildcards) else f"{input.fq1} {input.fq2}",
#         single = lambda wildcards: '' if isPE(wildcards) else '--single',
#         frag_length = lambda wildcards: '' if isPE(wildcards) else f"-l {se_frag_length}",
#         frag_sd = lambda wildcards: '' if isPE(wildcards) else f"-s {se_frag_sd}"
#     log:
#         "logs/kallisto/{sample_id}.log"
#     threads: qthreads
#     singularity: kallisto_env
#     shell:
#         """
#         kallisto quant -i {input.idx} -o {output} -t {threads} \
#             {params.strand} {params.single} {params.frag_length} {params.frag_sd} \
#             {params.fqs} &> >(tee {log})
#         """

rule salmon:
    input:
        unpack(get_fqs),
        idx=salmon_idx
    output:
        directory("salmon/{sample_id}")
    params:
        fqs = lambda wildcards, input: f"-r {input.fq}" if not isPE(wildcards) else f"-1 {input.fq1} -2 {input.fq2}",
        reps = lambda wildcards, input: f"--numGibbsSamples 50 -d" if salmon_add_reps else ""
    log:
        "logs/salmon/{sample_id}.log"
    threads: qthreads
    singularity: salmon_env 
    shell:
        """
        salmon quant -i {input.idx} -l A {params.fqs} -p {threads} -o {output} {params.reps}
        """

rule deseq2_init:
    input:
        unpack(get_quants),
        samples=samples_fp
    output:
        deseq="deseq2/all.rds",
        cts="results/normalized_counts.rds",
        se_tx="deseq2/se_tx.rds",
    params:
        aligner=quant_program,
        formula=design_formula,
        levels=var_levels
    log:
        "logs/deseq2/init.log"
    threads: dthreads
    singularity: r_env 
    script:
        "scripts/deseq2.R"

rule diffexp:
    input:
        gene="deseq2/all.rds"
    output:
        mleres = "results/{contrast}/mle_foldchanges.csv",
        mapres = "results/{contrast}/map_foldchanges.csv",
        mlema = "results/{contrast}/mle_ma.svg",
        mapma = "results/{contrast}/map_ma.svg",
    params:
        formula=design_formula,
        contrast = get_contrast
    log:
        "logs/deseq2/{contrast}.log"
    threads: dthreads
    singularity: r_env
    script:
        "scripts/diffexp.R"

rule multiqc:
    input:
        get_multiqc_input
    output:
        "qc/multiqc_report.html"
    log:
        "logs/multiqc.log"
    singularity: multiqc_env 
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
    singularity: r_env
    script:
        "scripts/plot_pca.R"

rule fastqc:
    input:
        unpack(get_fqs)
    output:
        directory("qc/fastqc/{sample_id}")
    singularity: fastqc_env
    log:
        "logs/fastqc/{sample_id}.log"
    shell:
        """
        tmpdir=qc/fastqc/.{wildcards.sample_id}.tmp
        mkdir $tmpdir
        mkdir {output}
        fastqc {input} -o {output} &> >(tee {log}) -d $tmpdir
        rm -r $tmpdir
        """


rule validateFastq:
    input:
        unpack(get_fqs)
    output:
        "qc/validateFastq/{sample_id}.txt"
    params:
        fqs = lambda wildcards, input: f"-i {input.fq}" if not isPE(wildcards) else f"-i {input.fq1} -j {input.fq2}"
    singularity: validateFastq_env
    shell:
        """
        biopet-validatefastq {params} 2>&1 | tee {output}
        """

rule checkFastqValidation:
    input:
        expand("qc/validateFastq/{sample_id}.txt", sample_id = samples['id'])
    output:
        "qc/validateFastq_summary.csv"
    script:
        "scripts/summarize_fastqValidation.py"
