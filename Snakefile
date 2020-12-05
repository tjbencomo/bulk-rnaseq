import sys
import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version

configfile: "config.yaml"

samples = pd.read_csv(config['samples'], dtype=str)
samples['id'] = samples['patient'] + '-' + samples['sample']
samples = samples.set_index(["id"], drop=False)
samples = samples.sort_index()

kallisto_idx = config['kallisto_index']
kallisto_tx2g = config['kallisto_tx2gene']
kallisto_threads = config['kallisto_threads']

INFER_READS = 50000

def get_fqs(wildcards):
    return {
        'fq1' : samples.loc[(wildcards.sample_id), 'fq1'],
        'fq2' : samples.loc[(wildcards.sample_id), 'fq2']
        }

def getStrand(wildcards):
    s = samples.loc[wildcards.sample_id, 'strandedness'].lower()
    if s == 'forward':
        return '--fr-stranded'
    elif s == 'reverse':
        return '--rf-stranded'
    elif s == 'none':
        return ''
    else:
        raise ValueError(f"Unrecognized strand type for {wildcards.sample_id}")

rule targets:
    input:
        expand("kallisto/{sample_id}", sample_id=samples['id'])

rule kallisto:
    input:
        unpack(get_fqs),
        idx=kallisto_idx
    output:
        directory("kallisto/{sample_id}")
    params:
        strand = lambda wildcards: getStrand(wildcards)
    log:
        "logs/kallisto/{sample_id}.log"
    threads: 1
    conda:
        "envs/quant.yml"
    shell:
        """
        kallisto quant -i {input.idx} -o {output} -t {threads} {params.strand} \
            {input.fq1} {input.fq2} | tee {log}
        """
