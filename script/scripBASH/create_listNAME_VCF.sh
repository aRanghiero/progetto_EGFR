#!/bin/bash

#08/09/21

##### con questo script voglio rinominare i file VCF togliendo il codice istologico del paziente e sostituendolo con un numero progressivo

#vedo dove mi trovo
pwd
clear
# attivo l'environment di conda con il programma BCFTOOLS 

eval "$(conda shell.bash hook)"
conda activate PicardTools

mkdir RenameVCF
mkdir NameFolder
PICARD="/home/albe/miniconda3/envs/PicardTools/share/picard-2.26.1-0/picard.jar" #set come variabile il percorso per utilizzare il file picard.jar

## creo un ciclo FOR per esportare i nomi dei campioni
# cd NameFolder
for file in *.vcf
do
	me=$(grep  -v "##" "$file" |  grep -v "chr" | awk -F ' ' '{print $10}') # prende l'id del paziente
	echo $me
	java -jar $PICARD MakeVcfSampleNameMap --INPUT "$file"  --OUTPUT NameFolder/"$(echo $file)".csv --VERBOSITY ERROR
done

###### questa parte la devo eseguire in bash
# cat *.csv > name_list.csv #unisco tutti i file csv
# devo creare la numerazione alternativa
# cat name_list.csv | tr  ',' '\t' > name_list.txt # converto il csv in txt
### dopo di che uso questo file per rinominare i file vcf


