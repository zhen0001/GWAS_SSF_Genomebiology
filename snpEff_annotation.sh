#!/bin/bash
#SBATCH -A SNIC2017-7-296
#SBATCH -p core -n 3 
#SBATCH -t 8:00:00
#SBATCH --mail-user zhiqiang.chen@slu.se
#SBATCH --mail-type=ALL

module load bioinfo-tools 
module load java
module load  htslib

folder="/proj/uppstore2017112/b2017108_nobackup/GenomicSelection/Viafull_GenoSel_Chen/GWASSSF"

# annotation of 5378 indivudals 
# 2019 -09 -07 
# after annotated file, which is not the gz file,  so we don't zat the file for next step

#  zcat $folder/SSF_indel_Bial_mis0.3.d140_GQ10_mis0.7.maf0.005_HWELocus.Imputed.sorted.P5238.vcf.gz |  java -Xmx10g -jar snpEff.jar  picea_abies  > $folder/SSF_indel_Bial_mis0.3.d140_GQ10_mis0.7.maf0.005_HWELocus.Imputed.sorted.P5238.ann.vcf


 java -Xmx10g -jar SnpSift.jar   extractFields $folder/SSF_indel_Bial_mis0.3.d140_GQ10_mis0.7.maf0.005_HWELocus.Imputed.sorted.P5238.ann.vcf  CHROM POS REF ALT "ANN[*].EFFECT"  > $folder/SSF_indel_Bial_mis0.3.d140_GQ10_mis0.7.maf0.005_HWELocus.Imputed.sorted.P5238.ann.txt
