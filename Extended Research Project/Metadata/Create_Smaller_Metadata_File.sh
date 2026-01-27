#This file has the code which was used to cut down the original larger CSV metadata file to have only the SRR_IDs and the ReadType
cut -d ',' -f 1, 26 SraRunTable.csv | tail -n +2 > SRR_ReadType.csv

#The column "LibraryLayout", which contains information about whether the reads were single or paired-end, was selected using column 26, despite being in column 19, because cut reads every comma in the file.  
