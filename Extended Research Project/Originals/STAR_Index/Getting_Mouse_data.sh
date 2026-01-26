#This is the code used to get the GteX and Fasta file that were used to generate the index
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38/gencode.vM38.primary_assembly.annotation.gtf.gz #for the GteX File
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38/GRCm39.primary_assembly.genome.fa.gz #For FastA file

#The link where the data was found was: https://www.gencodegenes.org/mouse/
#Then we need to use gunzip to be unzip both of the downloaded files before being about to use them for indexing

gunzip gencode.vM38.primary_assembly.annotation.gtf.gz
gunzip GRCm39.primary_assembly.genome.fa.gz
