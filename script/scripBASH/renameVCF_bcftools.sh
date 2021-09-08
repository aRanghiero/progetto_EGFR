#!/bin/bash

#08/09/21

##### con questo script voglio rinominare i file VCF togliendo il codice istologico del paziente e sostituendolo con un numero progressivo

#vedo dove mi trovo
pwd
clear
# attivo l'environment di conda con il programma BCFTOOLS 

eval "$(conda shell.bash hook)"
conda activate bcf_sam_vcf_tolls

## creo un ciclo FOR per cambiare i nomi dei campioni

for file in *.vcf
do
	me=$(grep  -v "##" "$file" |  grep -v "chr" | awk -F ' ' '{print $10}') # prende l'id del paziente
	echo $me
	bcftools reheader -s NameFolder/name_list.txt "$file" -o RenameVCF/"$(echo $file | sed 's/Non-Filtered/Rename_xxx/g')" --threads 5
	
done


