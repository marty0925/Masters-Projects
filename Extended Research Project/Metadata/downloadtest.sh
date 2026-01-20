#!/bin/bash -l
#SBATCH --job-name=sra_download
#SBATCH --output=logs/slurm_%A_%a.out
#SBATCH --time=48:00:00
#SBATCH --mem=16G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

# Note to future me: run this script using 'sbatch --array=1-500%20 downloadSRA.sh'
# And make sure to be in the same folder as the SRR_Id.txt file

# Load the SRA module
module load sra-tools/3.0.3-gcc-13.2.0

# grabs the SRR ID using the array task ID
SRR_ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" SRR_IDs.txt)
echo "Processing $SRR_ID (Task ID: $SLURM_ARRAY_TASK_ID)"

# Adding the Data Output Folder
OUTPUT_DIR="All_FastQ_Files"
mkdir -p "$OUTPUT_DIR"

# Prefetch to download the .SRA file temporarily in the current folder
prefetch --output-directory . "$SRR_ID"

# Check the files downloaded without any issues and convert the SRA file into a fastQ file
if [ $? -eq 0 ]; then
    echo "Validation passed. Converting..."
    
    fastq-dump --split-3 --outdir "$OUTPUT_DIR" "${SRR_ID}/${SRR_ID}.sra"
    
    if [ $? -eq 0 ]; then
        echo "Success. Cleaning up .sra file."
        
        rm -rf "$SRR_ID" #Remove the temporary folder containing the .sra file
    else
        echo "ERROR: Conversion failed for $SRR_ID" >&2
        exit 1
    fi
else
    echo "ERROR: Download failed for $SRR_ID" >&2
    exit 1
fi
