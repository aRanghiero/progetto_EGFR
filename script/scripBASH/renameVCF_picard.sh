#!/bin/bash

#07/09/21

##### con questo script voglio rinominare i file VCF togliendo il codice istologico del paziente e sostituendolo con un numero progressivo

#vedo dove mi trovo
pwd
clear
# attivo l'environment di conda con il programma Picard 

eval "$(conda shell.bash hook)"
conda activate PicardTools 

mkdir RenameVCF
mkdir RenameVCF/oldName # creo una cartella dove salvare i vcf rinominati

i=1
PICARD="/home/albe/miniconda3/envs/PicardTools/share/picard-2.26.1-0/picard.jar" #set come variabile il percorso per utilizzare il file picard.jar

# faccio eseguire un ciclo FOR per rinominare i file VCF

for file in *.vcf
do
	me=$(grep  -v "##" "$file" |  grep -v "chr" | awk -F ' ' '{print $10}') # prende l'id del paziente
	echo $me
	i=$((i+1)) #creo il numero progressivo da sostituire all'id del paziente
	echo $i
	java -jar $PICARD MakeVcfSampleNameMap --INPUT "$file"  --OUTPUT RenameVCF/oldName/"$(echo $me)".csv --VERBOSITY ERROR
	java -jar $PICARD RenameSampleInVcf -I "$file" -O RenameVCF/"$(echo $file | sed 's/Non-Filtered/Rename_xxx/g')" --NEW_SAMPLE_NAME $i --VERBOSITY ERROR
	
done

