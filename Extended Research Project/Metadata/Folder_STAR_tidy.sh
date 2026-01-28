#!/bin/bash
#SBATCH --job-name=Folder_Creation
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --time=08:00:00

SRR_ID="/scratch/prj/dn_hamidlab_msc_proj/k25115470/metadata/SRR_IDs.txt"
Location="/scratch/prj/dn_hamidlab_msc_proj/k25115470/Results/Human_Star_Alignment/Mapped_results"

while IFS=',' read -r SRR_ID_Col; do

	FinLoc="/scratch/prj/dn_hamidlab_msc_proj/k25115470/Results/Human_Star_Alignment/Mapped_results/$SRR_ID_Col"

	if [ -d "$FinLoc" ]; then
		echo "Directory already exists"
else
	echo "Making Folder for $SRR_ID"
		mkdir -p "$FinLoc"
		mv "$Location"/"$SRR_ID_Col"_* "$Location"/"$SRR_ID_Col"/
fi
done < "$SRR_ID"
