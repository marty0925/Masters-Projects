#!/bin/bash
#SBATCH --job-name=sra_parallel
#SBATCH --output=logs/parallel_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20    
#SBATCH --mem=64G 
#SBATCH --time=08:00:00

module load sra-tools/3.0.3-gcc-13.2.0
module load parallel/20220522-gcc-13.2.0

#Changing settings for parallel
export NCBI_SETTINGS="/tmp/ncbi_settings_${SLURM_JOB_ID}.mkfg"

# Using the GNU Parallel
cat SRR_IDs.txt |tr -d '\r' | parallel -j $SLURM_CPUS_PER_TASK --joblog download_log.txt  \
    "fastq-dump --split-files --gzip --outdir All_FastQ_Files {}" # --joblog to keep a log of which files succeeded/failed and outputting to the All_FastQ_Files Directory

rm -f "$NCBI_SETTINGS" #To remove the setting folder
