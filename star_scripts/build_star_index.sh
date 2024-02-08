#!/bin/bash

#SBATCH --cpus-per-task=24
#SBATCH --mem=128G
#SBATCH --job-name=star_index
#SBATCH --time=06:00:00
#SBATCH --output=star_%j.out
#SBATCH --error=star_%j.err

module load StdEnv/2020
module load gcc/9.3.0
module load star/2.7.9a


FASTA="/home/umaiyal1/scratch/refdata-gex-GRCh38-2020-A/fasta/genome.fa"
GTF="/home/umaiyal1/scratch/refdata-gex-GRCh38-2020-A/genes/genes.gtf"

STAR --runThreadN 24 --runMode genomeGenerate --genomeDir /home/umaiyal1/scratch/star_index --genomeFastaFiles $FASTA --sjdbGTFfile $GTF
