#!/bin/bash

## Script for sorting and indexing bam files
## Date: 12 Nov 2020 
##
## Example usage:
## inDir=/Shares/CL_Shared/data/atma/cohen2017_chip/bams sbatch --array 0-42 indexAndSortBam.q

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --time=1-00:00
#SBATCH --mem=64GB

# Job name and output
#SBATCH -J indexAndSortBam
#SBATCH -o /Users/%u/slurmOut/slurm-%A_%a.out
#SBATCH -e /Users/%u/slurmErr/slurm-%A_%a.err

# define query bam files
queries=($(ls $inDir/*.bam | xargs -n 1 basename))

# load modules
module load samtools

# run the thing
pwd; hostname; date

echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}
echo $(date +"[%b %d %H:%M:%S] Sorting bam file...")

samtools sort ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]} -o ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}.sorted

echo $(date +"[%b %d %H:%M:%S] Indexing bam file...")

samtools index ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}.sorted

echo $(date +"[%b %d %H:%M:%S] Done!")
