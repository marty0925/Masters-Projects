#!/bin/bash
#SBATCH --job-name=STAR_Index 
#SBATCH --ntasks=1  
#SBATCH --mem=16G 
#SBATCH --cpus-per-task=25   
#SBATCH --time=08:00:00

module load star/2.7.10b-gcc-13.2.0
STAR --runThreadN 25 \
--runMode genomeGenerate \
--genomeDir /scratch/prj/dn_hamidlab_msc_proj/k25115470/Originals/RefGenome_GRCm39/Index \
--genomeFastaFiles /scratch/prj/dn_hamidlab_msc_proj/k25115470/Originals/RefGenome_HumanGRCh38.p14/GRCh38.primary_assembly.genome.fa \
--sjdbGTFfile /scratch/prj/dn_hamidlab_msc_proj/k25115470/Originals/RefGenome_HumanGRCh38.p14/gencode.v49.primary_assembly.annotation.gtf  \
--sjdbOverhang 51

#At the end of the Log.out file, there will be a message similar to "Jan 26 11:47:11 ..... finished successfully"
