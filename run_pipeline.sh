#!/bin/bash
#SBATCH --job-name=
#SBATCH --output=
#SBATCH --nodes=1
#SBATCH --time=00-12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=6000
#SBATCH --mail-type=END
#SBATCH --mail-user=

set -e
cd $(pwd)
echo "Starting up snakemake..."
snakemake --cluster-config cluster.json -j 499 \
    --use-singularity \
    --cluster 'sbatch -p {cluster.partition} -t {cluster.time} --mem {cluster.mem} -c {cluster.ncpus} -o {cluster.out}'
