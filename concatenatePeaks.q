#!/bin/bash

## Script for merging peaks from all samples in a batch
## Date: 12 Nov 2020 
##
## Example usage:
## inDir=/Users/CL_Shared/data/atma/cohen2017_chip/filtered_peaks projectName=cohen2017_chip sbatch mergePeaks.q

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=1:00:00
#SBATCH --mem=8GB

# Job name and output
#SBATCH -J mergePeaks
#SBATCH -o /Users/%u/slurmOut/slurm-%j.out
#SBATCH -e /Users/%u/slurmErr/slurm-%j.err

# load modules
module load bedtools

# run the thing
pwd; hostname; date

echo "Merging peaks from all samples..."

# go to dir
cd $inDir

# concatenate all broadpeak regions from all samples
cat *.broadPeak.filtered > "$projectName"_all.tmp

# sort and merge intervals within 100bp from each other
bedtools sort -i "$projectName"_all.tmp > "$projectName"_all_sorted.tmp
bedtools merge -i "$projectName"_all_sorted.tmp -d 100 > "$projectName"_all_sorted_merged100bp.tmp

# label the regions
cat "$projectName"_all_sorted_merged100bp.tmp \
| awk '{print $0 "\t" "region"NR}' \
> ALL_"$projectName".broadPeaks

# remove tmp files
rm *.tmp

echo $(date +"[%b %d %H:%M:%S] Done!")
