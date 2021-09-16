#!/bin/bash

# 14/09/21

eval "$(conda shell.bash hook)"
conda activate bcf_sam_vcf_tolls 

#vedo dove mi trovo
pwd
cd RenemaVCF
mkdir ../cnv # creo una cartella dove le informazioni legate ai CNV in formato .csv

# faccio eseguire un ciclo FOR per estrarre i CNV dal file VCF
for file in *.vcf
do
	me=$(grep  -v "##" "$file" |  grep -v "chr" | awk -F ' ' '{print $10}') # prende l'id del paziente
	echo $me
	bcftools query -f '%CHROM\t%POS\t%ID\t%ALT\t[%CN\t]\t[%CDF_MAPD\t]\n' "$file" |awk '$5>=4'| sed 's/,/\t/g' | sed 's/:/\t/g'  | awk -v me="$me" '$15>=4 {print $3 "," me ","   "Amp"}' > cnv/"$(echo ${file%.vcf}.csv | sed 's/xxx/CNV/g')"
done
#cd cnv
#cat *.csv > tot_cnv.csv

