#Download the "SraRunRable.csv" file from the SRA Run Selector website
#Project ID: PRJNA940433 / GSE226482
#SRA Link: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA940433&o=acc_s%3Aa

 cut -d ',' -f 1 SraRunTable.csv | tail -n +2 > SRR_IDs.txt #Using the cut and tail commands to isolate the SRR IDs and generate a .txt file with all 1415 IDs. 
