#This file contains the code used to run and validate the dry run of the STAR Alignment conducted. 
#!/bin/bash
#SBATCH --job-name=STAR_Index
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --time=08:00:00

module load star/2.7.10b-gcc-13.2.0

GENOME_DIR="/scratch/prj/dn_hamidlab_msc_proj/k25115470/Originals/RefGenome_HumanGRCh38.p14/Index"
INPUT_CSV="/scratch/prj/dn_hamidlab_msc_proj/k25115470/metadata/SRR_ReadType.csv"
OUTPUT_DIR="/scratch/prj/dn_hamidlab_msc_proj/k25115470/Results/Human_Star_Alignment/Mapped_results"
FASTQ_DIR="/scratch/prj/dn_hamidlab_msc_proj/k25115470/Originals/All_FastQ_Files"
THREADS=16

mkdir -p "$OUTPUT_DIR"

# Loop through the CSV file
while IFS=',' read -r SRR_ID ReadType; do # read -r reads Run_id and Layout into variables

    echo "Processing sample: $SRR_ID"
    echo "Type: $ReadType"

    # Defining the output prefix
    OUT_PREFIX="${OUTPUT_DIR}/${SRR_ID}_"

    if [ "$ReadType" == "PAIRED" ]; then #Paried-end Side

        R1="${FASTQ_DIR}/${SRR_ID}_1.fastq.gz"
        R2="${FASTQ_DIR}/${SRR_ID}_2.fastq.gz"

        echo "Files: $R1 and $R2"

        echo "STAR --runThreadN "$THREADS" \
             --genomeDir "$GENOME_DIR" \
             --readFilesIn "$R1" "$R2" \
             --readFilesCommand zcat \
             --outFileNamePrefix "$OUT_PREFIX" \
             --outSAMtype BAM SortedByCoordinate \
             --quantMode GeneCounts"

    elif [ "$ReadType" == "SINGLE" ]; then #Single side
        R1="${FASTQ_DIR}/${SRR_ID}_1.fastq.gz"

        echo "File: $R1"

        echo "STAR --runThreadN "$THREADS" \
             --genomeDir "$GENOME_DIR" \
             --readFilesIn "$R1" \
             --readFilesCommand zcat \
             --outFileNamePrefix "$OUT_PREFIX" \
             --outSAMtype BAM SortedByCoordinate \
             --quantMode GeneCounts"

    else
        echo "ERROR: The Read Type '$ReadType' was not recognized for $SRR_ID. It will be Skipped."
    fi

done < "$INPUT_CSV"

#To validate the results, I used the grep command on the slurm.out file
grep "PAIRED" DRYRUN_slurm-31390695.out | wc # 616, 1232, 8008
grep "SINGLE" DRYRUN_slurm-31390695.out | wc # 799, 1598, 10387
# 616+799=1415, which is the same number of reads that we had. 





