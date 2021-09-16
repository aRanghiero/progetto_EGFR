#!/bin/bash

## cerchiamo di creare uno script per convertire i VCF in MAf
### funziona il 11/05/21
# il programma si trova https://github.com/mskcc/vcf2maf


#vedo dove mi trovo
pwd

# attivo l'environment di conda con il programma vcf2maf 
# conda activate VCF2MAF 
eval "$(conda shell.bash hook)"
conda activate VCF2MAF 
# mi devo trovare dentro la cartella del porgramma /media/albe/Volume/master/tirocinio/mskcc-vcf2maf-754d68a

# creo delle cartelle
cd vcf_tot # dentro questa cartella si devono trovare i file VCF
mkdir maf # creo una cartella dove salvare i maf
mkdir maf/temp # creo una cartella temporanea dove vengono salvati i file temp

# faccio eseguire un ciclo FOR per convertire tutti i VCF in MAF
for file in *.vcf
#           ^  Bash esegue l’espansione del nome del file
#+             nelle espressioni che riconoscono il globbing.

do
	me=$(grep  -v "##" "$file" |  grep -v "chr" | awk -F ' ' '{print $10}')
	echo "$file"
	echo $me
	perl ../vcf2maf.pl --input-vcf "$file" --output maf/"${file%.vcf}.maf"  --vep-path /home/albe/miniconda3/envs/VCF2MAF/bin  --vep-overwrite --tumor-id $me --ref-fasta ../hg19.fa --retain-fmt AO,AF --tmp-dir maf/temp 
	

done # chiude il ciclo for





