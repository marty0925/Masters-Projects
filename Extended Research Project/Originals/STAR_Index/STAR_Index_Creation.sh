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
--genomeFastaFiles /scratch/prj/dn_hamidlab_msc_proj/k25115470/Originals/RefGenome_GRCm39/GRCm39.primary_assembly.genome.fa \
--sjdbGTFfile /scratch/prj/dn_hamidlab_msc_proj/k25115470/Originals/RefGenome_GRCm39/gencode.vM38.primary_assembly.annotation.gtf \
--sjdbOverhang 51
