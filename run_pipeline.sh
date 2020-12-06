#!/bin/bash
#SBATCH --job-name=bulk-rnaseq
#SBATCH --output=
#SBATCH --nodes=1
#SBATCH --time=01-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4000
#SBATCH --mail-type=END
#SBATCH --mail-user=

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}"  )" >/dev/null 2>&1 && pwd  )"
set -e
cd $DIR
snakemake --cluster-config cluster.json -j 499 \
    --use-conda --use-singularity \
    --cluster 'sbatch -p {cluster.partition} -t {cluster.time} --mem {cluster.mem} -c {cluster.ncpus} -o {cluster.out}'
