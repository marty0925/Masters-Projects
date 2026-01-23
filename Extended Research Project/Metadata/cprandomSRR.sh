#!/bin/bash
#SBATCH --job-name=FastQc_20_choosen
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20    
#SBATCH --mem=16G 
#SBATCH --time=02:00:00

#To define the paths for later
ID_LIST="Random_SRR_IDs.txt"
SOURCE_DIR="/scratch/prj/dn_hamidlab_msc_proj/k25115470/Originals/All_FastQ_Files"
DEST_DIR="/scratch/prj/dn_hamidlab_msc_proj/k25115470/Originals/Fastqc"

mkdir -p "$DEST_DIR"

#the Loop itself
cat "$ID_LIST" | tr -d '\r' | while read -r SRR_ID; do

    # Skip empty lines
    if [ -z "$SRR_ID" ]; then continue; fi
    
    FILES_TO_COPY=$(ls "$SOURCE_DIR/"*"$SRR_ID"* 2>/dev/null) #Use wildcard to select the SRR IDs with extension

    if [ -n "$FILES_TO_COPY" ]; then
        echo "Found matches for $SRR_ID. Copying..."
        
        #copy any files which match that pattern
        cp "$SOURCE_DIR/"*"$SRR_ID"* "$DEST_DIR/"
        
    else
        echo "Error: No files found containing ID: $SRR_ID"
    fi

done
