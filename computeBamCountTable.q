#!/bin/bash

## Script for computing bam read count table
## Date: 12 Nov 2020 
##
## Example usage:
## bamOrder=/Shares/CL_Shared/data/atma/cohen2017_chip/mybams.txt 
## bed=/Shares/CL_Shared/data/atma/cohen2017_chip/ALL_peaks.narrowPeak 
## outDir=/Shares/CL_Shared/data/atma/cohen2017_chip/bamCounts 
## sbatch computeBamCountTable.q
## 
## Note:
## bamOrder is a text file containing the full path and file names of all the bams, in the order they should appear as columns in the table
## You need to create the bamOrder file before running this script
## And you must include the full path to the bamOrder file, as shown above
## 
## Here's what the inside of a bamOrder file should look like:
## /Shares/CL_Shared/data/atma/cohen2017_chip/file1.bam.sorted
## /Shares/CL_Shared/data/atma/cohen2017_chip/file2.bam.sorted
## /Shares/CL_Shared/data/atma/cohen2017_chip/file3.bam.sorted
## /Shares/CL_Shared/data/atma/cohen2017_chip/file4.bam.sorted
## 

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=1-00:00
#SBATCH --mem=64GB

# Job name and output
#SBATCH -J bamCountTable
#SBATCH -o /Users/%u/slurmOut/slurm-%j.out
#SBATCH -e /Users/%u/slurmErr/slurm-%j.err

# load modules
module load bedtools

# define key variables
bams=$(cat $bamOrder)

# run the thing
pwd; hostname; date

echo "Starting bedtools..."
echo $(date +"[%b %d %H:%M:%S] Computing matrix of bam read counts...")

bedtools multicov -bams ${bams} -bed ${bed} > $outDir/bamCountsWithinPeaks.tab

echo $(date +"[%b %d %H:%M:%S] Done!")
