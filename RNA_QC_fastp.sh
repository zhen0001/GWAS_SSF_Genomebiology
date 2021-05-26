#!/bin/bash
#SBATCH -A SNIC2017-7-296
#SBATCH -p core -n 4
#SBATCH -t 40:00:00
#SBATCH --mail-user zhiqiang.chen@slu.se
#SBATCH --mail-type=ALL

module load bioinfo-tools
module load vcftools
module load vcflib/2017-04-04

# Un this study, I used fastp to run quality control and provide clean data

# just example for single sample
#  fastp -i ../RNA/L9_5_1.fq.gz  -o ../Fastp_outdir/L9_5_1.fq.gz  -I  ../RNA/L9_5_2.fq.gz  -O ../Fastp_outdir/L9_5_2.fq.gz    -j $i.json -h $j.html



# using loop for fastq 

for i in $(ls ../RNA/*1.fq.gz)
do 
name=${i%_*}
fastp -i ${i%_*}_1.fq.gz  -o ../Fastp_outdir/${name##*/}_1.fq.gz  -I  ${i%_*}_2.fq.gz  -O ../Fastp_outdir/${name##*/}_2.fq.gz     -j ../Fastp_outdir/${name##*/}.json -h ../Fastp_outdir/${name##*/}.html

done


