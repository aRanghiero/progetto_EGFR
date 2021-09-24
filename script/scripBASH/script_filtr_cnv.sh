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

# cat 3_Rename_CNV_2019-12-10_23_48_12.csv 4_RNA_v2_Rename_CNV_2019-12-10_23_50_06.csv 9_RNA_v1_Rename_CNV_2019-12-10_23_54_52.csv 10_RNA_v1_Rename_CNV_2019-12-10_23_57_51.csv 15_RNA_v1_Rename_CNV_2019-12-11_00_02_47.csv 18_RNA_v1_Rename_CNV_2019-12-11_00_09_51.csv 20_RNA_v1_Rename_CNV_2019-12-11_00_13_11.csv 21_RNA_v1_Rename_CNV_2019-12-11_00_13_40.csv 22_RNA_v1_Rename_CNV_2019-12-11_00_25_34.csv 23_RNA_v1_Rename_CNV_2019-12-11_00_25_59.csv 25_RNA_v1_Rename_CNV_2019-12-11_00_28_53.csv 32_RNA_v1_Rename_CNV_2019-12-11_00_37_24.csv 34_RNA_v1_Rename_CNV_2020-01-20_01_38_07.csv 41_RNA_v1_Rename_CNV_2020-01-31_01_26_25.csv 42_RNA_v1_Rename_CNV_2020-02-03_04_43_13.csv 46_RNA_v1_Rename_CNV_2020-02-07_06_54_49.csv 50_Rename_CNV_2020-02-24_05_56_59.csv 51_Rename_CNV_2020-08-05_01_39_35.csv 52_Rename_CNV_2020-08-05_01_41_05.csv 53_Rename_CNV_2020-08-05_01_01_27.csv 54_Rename_CNV_2020-08-05_01_14_16.csv 56_Rename_CNV_2020-08-05_01_14_59.csv 57_Rename_CNV_2020-08-05_01_21_36.csv 58_Rename_CNV_2020-08-05_01_21_58.csv 59_Rename_CNV_2020-08-05_01_34_03.csv 60_Rename_CNV_2020-08-05_01_22_41.csv 61_Rename_CNV_2020-08-06_00_28_37.csv 62_Rename_CNV_2020-08-06_00_28_54.csv 63_Rename_CNV_2020-08-06_00_29_16.csv 65_Rename_CNV_2020-08-06_00_30_44.csv 66_Rename_CNV_2020-08-06_00_31_28.csv 67_Rename_CNV_2020-08-06_00_31_49.csv 68_Rename_CNV_2020-08-06_00_32_24.csv 20-C-002683_v1_20-C-002683_RNA_v1_Rename_CNV_2021-06-28_04_10_20.csv 20-M-001548_v1_20-M-001548_RNA_v1_Rename_CNV_2021-06-28_04_17_30.csv 20-I-013287b_v1_20-I-013287b_RNA_v1_Rename_CNV_2021-06-28_04_35_44.csv 20-C-002911_v1_20-C-002911_RNA_v1_Rename_CNV_2021-06-28_04_36_02.csv 20-C-003624_v1_20-C-003624_RNA_v1_Rename_CNV_2021-06-28_04_37_14.csv 20-M-002416_v1_20-M-002416_RNA_v1_Rename_CNV_2021-06-28_04_42_44.csv 20-I-018638_v1_20-I-018638_RNA_v1_Rename_CNV_2021-06-28_04_50_12.csv 20-I-019632_v1_20-I-019632_RNA_v1_Rename_CNV_2021-06-28_04_50_56.csv 21-I-000209_v1_21-I-000209_RNA_v1_Rename_CNV_2021-06-28_04_51_18.csv 20-I-021069_v1_20-I-021069_RNA_v1_Rename_CNV_2021-06-28_04_51_49.csv 21-I-000031_v1_21-I-000031_RNA_v1_Rename_CNV_2021-06-28_04_58_12.csv 21-C-000276_v1_21-C-000276_RNA_v1_Rename_CNV_2021-06-28_04_59_14.csv 21-I-000404_v1_21-I-000404_RNA_v1_Rename_CNV_2021-06-28_05_04_23.csv 21-C-000569_v1_21-C-000569_RNA_v1_Rename_CNV_2021-06-28_05_04_23.csv 15-I-000801_v1_15-I-000801_RNA_v1_Rename_CNV_2021-06-28_05_07_37.csv 21-I-003450_v1_21-I-003450_RNA_v1_Rename_CNV_2021-06-28_05_07_42.csv 21-C-001140_v1_21-C-001140_RNA_v1_Rename_CNV_2021-06-28_05_11_15.csv 21-C-001672_v1_21-C-001672_RNA_v1_Rename_CNV_2021-06-28_05_12_31.csv 21-C-002073_v1_21-C-002073_RNA_v1_Rename_CNV_2021-06-29_01_31_42.csv 21-I-009552_v2_21-I-009552_RNA_v2_Rename_CNV_2021-06-29_01_32_01.csv 21-C-001917_v1_21-C-001917_RNA_v1_Rename_CNV_2021-06-29_01_36_33.csv 21-I-009074_v1_21-I-009074_RNA_v1_Rename_CNV_2021-06-29_01_36_55.csv 21-C-002375_v1_21-C-002375_RNA_v1_Rename_CNV_2021-06-29_01_42_20.csv > CNV_EX21.csv

# cat 2_RNA_v1_Rename_CNV_2019-12-10_23_39_36.csv 5_RNA_v1_Rename_CNV_2019-12-10_23_50_35.csv 6_RNA_v1_Rename_CNV_2019-12-10_23_50_42.csv 7_RNA_v1_Rename_CNV_2019-12-10_23_54_16.csv 8_RNA_v1_Rename_CNV_2019-12-10_23_54_42.csv 11_RNA_v1_Rename_CNV_2019-12-10_23_58_11.csv 12_RNA_v1_Rename_CNV_2019-12-10_23_58_44.csv 13_RNA_v1_Rename_CNV_2019-12-11_00_01_18.csv 14_RNA_v1_Rename_CNV_2019-12-11_00_02_17.csv 16_RNA_v1_Rename_CNV_2019-12-11_00_09_14.csv 17_RNA_v1_Rename_CNV_2019-12-11_00_09_33.csv 19_RNA_v1_Rename_CNV_2019-12-11_00_12_52.csv 24_RNA_v1_Rename_CNV_2019-12-11_00_26_10.csv 26_RNA_v1_Rename_CNV_2019-12-11_00_29_12.csv 27_RNA_v1_Rename_CNV_2019-12-11_00_30_21.csv 28_RNA_v1_Rename_CNV_2019-12-11_00_34_32.csv 29_RNA_v1_Rename_CNV_2019-12-11_00_34_34.csv 30_RNA_v1_Rename_CNV_2019-12-11_00_35_05.csv 31_RNA_v1_Rename_CNV_2019-12-11_00_37_08.csv 33_RNA_v1_Rename_CNV_2019-12-11_00_37_26.csv 35_RNA_v2_Rename_CNV_2020-01-20_01_42_50.csv 36_RNA_v2_Rename_CNV_2020-01-20_01_51_25.csv 37_RNA_v1_Rename_CNV_2020-01-20_01_55_27.csv 38_RNA_v1_Rename_CNV_2020-01-31_01_06_46.csv 39_RNA_v1_Rename_CNV_2020-01-31_01_15_24.csv 40_RNA_v1_Rename_CNV_2020-01-31_01_24_58.csv 43_RNA_v1_Rename_CNV_2020-02-03_05_50_51.csv 44_RNA_v1_Rename_CNV_2020-02-07_06_50_53.csv 45_RNA_v1_Rename_CNV_2020-02-07_06_54_25.csv 47_RNA_v1_Rename_CNV_2020-02-14_03_18_02.csv 48_RNA_v1_Rename_CNV_2020-02-18_05_29_51.csv 49_RNA_v1_Rename_CNV_2020-02-24_05_56_34.csv 69_Rename_CNV_2020-08-06_01_08_16.csv 70_Rename_CNV_2020-08-06_01_08_34.csv 71_Rename_CNV_2020-08-06_01_08_59.csv 20-M-001487_v1_20-M-001487_RNA_v1_Rename_CNV_2021-06-28_04_17_54.csv 20-I-010810_v1_20-I-010810_RNA_v1_Rename_CNV_2021-06-28_04_18_30.csv 20-C-002372_v1_20-C-002372_RNA_v1_Rename_CNV_2021-06-28_04_18_35.csv 20-I-012720_v1_20-I-012720_RNA_v1_Rename_CNV_2021-06-28_04_36_49.csv 20-C-003274b_v1_20-C-003274b_RNA_v1_Rename_CNV_2021-06-28_04_37_00.csv 20-C-003994_v2_20-C-003994_RNA_v2_Rename_CNV_2021-06-28_04_42_03.csv 20-I-015201_v1_20-I-015201_RNA_v1_Rename_CNV_2021-06-28_04_41_57.csv 20-I-017812_v1_20-I-017812_RNA_v1_Rename_CNV_2021-06-28_04_43_02.csv 20-I-017106_v1_20-I-017106_RNA_v1_Rename_CNV_2021-06-28_04_47_12.csv 20-I-019355_v1_20-I-019355_RNA_v1_Rename_CNV_2021-06-28_04_51_15.csv 21-I-001080_v1_21-I-001080_RNA_v1_Rename_CNV_2021-06-28_04_58_07.csv 21-I-000899_v1_21-I-000899_RNA_v1_Rename_CNV_2021-06-28_04_58_28.csv 21-I-000644_v1_21-I-000644_RNA_v1_Rename_CNV_2021-06-28_04_58_34.csv 21-S-000015b_v1_Rename_CNV_2021-06-28_05_03_19.csv 20-C-005060_v1_20-C-005060_RNA_v1_Rename_CNV_2021-06-28_05_03_13.csv 21-C-000475_v1_21-C-000475_RNA_v1_Rename_CNV_2021-06-28_05_04_05.csv 21-C-000745_v1_21-C-000745_RNA_v1_Rename_CNV_2021-06-28_05_07_18.csv 19-I-017983_v1_19-I-017983_RNA_v1_Rename_CNV_2021-06-28_05_07_54.csv 21-C-001187_v1_21-C-001187_RNA_v1_Rename_CNV_2021-06-28_05_07_38.csv 21-C-000036_v1_21-C-000036_RNA_v1_Rename_CNV_2021-06-28_05_11_56.csv 21-C-001481_v1_21-C-001481_RNA_v1_Rename_CNV_2021-06-28_05_11_53.csv 21-I-006272_v1_21-I-006272_RNA_v1_Rename_CNV_2021-06-28_05_12_09.csv 21-I-007888_v1_21-I-007888_RNA_v1_Rename_CNV_2021-06-29_01_30_35.csv 21-C-001826_v1_21-C-001826_RNA_v1_Rename_CNV_2021-06-29_01_30_59.csv 21-C-001916_v1_21-C-001916_RNA_v1_Rename_CNV_2021-06-29_01_31_20.csv 21-M-001157_v1_21-M-001157_RNA_v1_Rename_CNV_2021-06-29_01_37_10.csv 21-I-008681_v1_21-I-008681_RNA_v1_Rename_CNV_2021-06-29_01_41_45.csv 21-C-002271_v2_21-C-002271_RNA_v2_Rename_CNV_2021-06-29_01_42_01.csv 21-I-010197_v1_21-I-010197_RNA_v1_Rename_CNV_2021-06-29_01_42_40.csv > CNV_EX19.csv