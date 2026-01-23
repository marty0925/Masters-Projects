#!/bin/bash
#SBATCH --job-name=FastQc_20_choosen
#SBATCH --ntasks=1    
#SBATCH --mem=16G 
#SBATCH --time=02:00:00

module load fastqc/0.12.1-gcc-13.2.0

fastqc /scratch/prj/dn_hamidlab_msc_proj/k25115470/Originals/Fastqc/*.gz  --outdir /scratch/prj/dn_hamidlab_msc_proj/k25115470/Originals/Fastqc/results
