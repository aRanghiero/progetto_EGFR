#!/bin/bash

#02/09/21

#vedo dove mi trovo
pwd

# attivo l'environment di conda con il programma bcftools 
# conda activate VCF2MAF 
eval "$(conda shell.bash hook)"
conda activate bcf_sam_vcf_tolls 

mkdir filt5% # creo una cartella dove salvare i vcf filtrati

# faccio eseguire un ciclo FOR per convertire tutti i VCF in MAF
for file in *.vcf
do
	me=$(grep  -v "##" "$file" |  grep -v "chr" | awk -F ' ' '{print $10}') # prende l'id del paziente
	echo $me
	# filtro al 5%
	bcftools view --min-af 0.050:nref "$file" --output-type v --output-file filt5%/"$(echo $file | sed 's/xxx/Filtered5/g')" --threads 5 
	

done
