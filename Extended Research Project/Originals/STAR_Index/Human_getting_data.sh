#This is the code used to get the GteX and Fasta file that were used to generate the human index
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.primary_assembly.annotation.gtf.gz #for the GteX File
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh38.primary_assembly.genome.fa.gz #For FastA file

#The link where the data was found was: https://www.gencodegenes.org/human/
#Then we need to use gunzip to be unzip both of the downloaded files before being about to use them for indexing

gunzip gencode.v49.primary_assembly.annotation.gtf.gz
gunzip GRCh38.primary_assembly.genome.fa.gz
