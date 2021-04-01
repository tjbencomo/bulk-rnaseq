#!/bin/bash
#SBATCH --job-name=bulk-rnaseq
#SBATCH --output=/scratch/users/tbencomo/logs/bulk-rnaseq.out
#SBATCH --nodes=1
#SBATCH --time=00-04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8000
#SBATCH --mail-type=END
#SBATCH --mail-user=tbencomo@stanford.edu

set -e
cd $1
echo "Starting up snakemake..."
snakemake --cluster-config cluster.json -j 499 \
    --use-conda --use-singularity \
    --cluster 'sbatch -p {cluster.partition} -t {cluster.time} --mem {cluster.mem} -c {cluster.ncpus} -o {cluster.out}'
