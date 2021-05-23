# Remove everything in the working environment
rm(list=ls())  

# only for the first time when install new R version
# BiocManager::install("http://www.bioconductor.org/biocLite.R")
# devtools::install_github("YaoZhou89/BLINK")
#  BiocManager::install("biocLite")

# if this is done, next time, we don't need it.  
# install.packages("gplots")
# install.packages("LDheatmap")
# install.packages("genetics")
# install.packages("ape")
# install.packages("EMMREML")
# install.packages("scatterplot3d")

# install.packages("bigmemory ")
# install.packages("biganalytics")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("multtest")) # this is only way to install package from BiocManager

library(bigmemory)
library(biganalytics)
library(multtest)
library(gplots)
library(LDheatmap)
library(genetics)
library(ape)
library(EMMREML)
library(compiler) #this library is already installed in R
library("scatterplot3d")

source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("/Users/zhen0001/Box\ Sync/Chen.function.R")

# The EMMA library was developed by Kang et al. (2008). One line was added to the library to handle the
# grid search with "NA" likelihood. The modified library can be installed by typing this command line:

source("http://zzlab.net/GAPIT/emma.txt")

setwd("/Users/zhen0001/BLINK-5238/")

#******************************************
# my own sample from southern Sweden
# C:/Users/zhen0001/Desktop/SSF8000/southern_population
#*******************************************
# rm(list=ls())   
require(data.table) # fread is similar to read.table, but faster and convenient all controls suchas sep, colcloasses, nrows,
data<-fread("out.012", na.strings=c("",".", "NA")) # "-1" it is great helpful for # mdata1[mdata1=="-1"]=NA
data[1:10, 1:10]; dim(data)
# data<-data[, 1:10001] # to get subset of data
data<-data[,-1]


indv<-fread("out.012.indv", header = FALSE) # the name of individuals
sortmarker_id<- fread("out.012.pos", header = FALSE) # the name of markers
length(unique(sortmarker_id$V1)) # to count the number of contig

colnames(data)<-paste(sortmarker_id$V1, sortmarker_id$V2, sep="_") # [1:10000]
row.names(data)<-indv$V1
myGD<-data.frame(indv$V1, data)
rm(data)
colnames(myGD)[1]="Taxa"
myGD[1:10, 1:10]

#write.table(myGD, file = "myData.dat", quote = FALSE, row.names = FALSE, col.names = FALSE, sep="\t")
# rm(myGD)
head(sortmarker_id)
# n<-dim(sortmarker_id[1:10000,])[1]
n<-dim(sortmarker_id)[1]

# SNP<-paste(sortmarker_id$V1, sortmarker_id$V2, sep="_")[1:10000]
SNP<-paste(sortmarker_id$V1, sortmarker_id$V2, sep="_")
length(SNP)
# myGM<-data.frame(SNP,rep(1:10, each=100), sortmarker_id$V2[1:1000])
myGM<-data.frame(SNP,c(rep(1:10, each=29000),rep(10, 9240)) , 1:n)
colnames(myGM)<-c("rs", "chr", "pos") # rs: reference snp
head(myGM)

# write.table(myGM, file = "/Users/zhen0001/BLINK-5238/myData.map", quote = FALSE, row.names = FALSE,col.names = TRUE, sep="\t" )
# rm(myGM)

#*****************************************************
# BLINK in R version
# association analysis
source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("http://zzlab.net/FarmCPU/FarmCPU_functions.txt")
library(BLINK)

# rm(x); rm(yc)
# genotype information data
# how to prepare this map 
#myGM=read.table("/Users/zhen0001/BLINK-5238/myData.map",head=T, col.names=c("rs", "chr", "pos") )
head(myGM); dim(myGM)
# myGM$chr=1

# myY_orignal <- read.table("/Users/zhen0001/Box\ Sync/Pheno_Geno_EBV_SSF_budset3_disease_Northphenology_Sävar_Frost_FFFF.txt", header = TRUE, sep="\t") # 2019_08_21
# head(myY_orignal)
# ROM<-droplevels(myY_orignal[myY_orignal$Group_FFF=="ROM" & !is.na(myY_orignal$Ajd_GV_MOE),])
# dim(ROM)
# myGD<-droplevels(myGD[is.element(myGD$Taxa,ROM$RAPiD_Genomics_Sample_Code), ])
# dim(myGD)
# rm(ROM)
# rm(myY_orignal)
# # genotype data
# myGD=read.big.matrix("/Users/zhen0001/BLINK-5238/myData.dat",header=F,sep="\t",type="char") 
# myGD<-myGD[, -1] # the first column should not appear in BLINK, but maybe in Gapit
myGD[1:10, 1:10]; dim(myGD)
myGD[1:]
# to get the frequency of each SNPs
freq<-function(x){
  (sum(x==2)*2+ sum(x==1))/(length(x)*2)
}

freq.snp<-apply(myGD[, -1],2, freq)
head(freq.snp[order(freq.snp)])
# hist(freq.snp)

str(freq.snp>=0.03)
sum(freq.snp>=0.03) # 135902
sum(freq.snp>=0.03 & freq.snp<=0.97 ) # 134605 # some snps is not the 2 is minor allele frequency
sum(freq.snp>=0.05 & freq.snp<=0.95 ) # 98187 # so there is something different

# to check if the MyGM$ rs == anno.snp$rs 
anno.snp<-read.table("/Users/zhen0001/BLINK-5238/SSF_indel_Bial_mis0.3.d140_GQ10_mis0.7.maf0.005_HWELocus.Imputed.sorted.P5238.ann.txt", header = TRUE, sep="\t", stringsAsFactors = T)
head(anno.snp); dim(anno.snp)
names(anno.snp)[5:30]<-paste("ANNEffect", 1:26, sep="_")
anno.snp$rs<-noquote(paste(anno.snp$CHROM, anno.snp$POS, sep="_"))
head(anno.snp); summary(anno.snp$rs)

# this should be correct for us
sum(is.element(anno.snp$rs, myGM$rs ))
identical(myGM$rs, as.factor(anno.snp$rs)) # whether it is identical


# it dependents what level I should do it. 
myGD<-myGD[, c(TRUE, freq.snp>=0.03 & freq.snp<=0.97)] # here I used 
myGM<-myGM[c(freq.snp>=0.03 & freq.snp<=0.97), ]

# myGD<-myGD[, c(freq.snp>=0.05 & freq.snp<=0.95)] # here I used 
# myGM<-myGM[c(freq.snp>=0.05 & freq.snp<=0.95), ]

# to get annotation for each snps
sub.anno.snp<-anno.snp[c(freq.snp>=0.03 & freq.snp<=0.97), ]

#*********************************
# to get all frequency, we may don't need it 
allele<-read.table("/Users/zhen0001/BLINK-5238/GWAS5238.frq.count", header = TRUE, sep = "\t")
head(allele); dim(allele)
allele$rs<-paste(allele$CHROM, allele$POS, sep="_")
sub.allele<-allele[c(freq.snp>=0.03 & freq.snp<=0.97), ]

# new freq.snp 133141
freq.snp<-apply(myGD[,-1],2, freq)
head(freq.snp); length(freq.snp)

# 
# mdat<-read.xls("/Users/zhen0001/Box\ Sync/GWAS/results/BLINK_produced_results/Table2/Table2_SNP_markers_F.xlsx", sheet = 3, header=T)
# head(mdat)

# the whole SSF phenotype data
# myY_orignal <- read.table("/Users/zhen0001/Box\ Sync/Pheno_Geno_EBV_SSF_budset3_disease_Northphenology.txt", header = TRUE, sep="\t") # 2019_08_21
myY_orignal <- read.table("/Users/zhen0001/Box\ Sync/Pheno_Geno_EBV_SSF_budset3_disease_Northphenology_Sävar_Ekebo11_F.txt", header = TRUE, sep="\t") # 2019_08_21
myY_orignal <- read.table("/Users/zhen0001/Box\ Sync/Pheno_Geno_EBV_SSF_budset3_disease_Northphenology_Sävar_Frost_F.txt", header = TRUE, sep="\t") # 2019_08_21
myY_orignal <- read.table("/Users/zhen0001/Box\ Sync/Pheno_Geno_EBV_SSF_budset3_disease_Northphenology_Sävar_Frost_FF.txt", header = TRUE, sep="\t") # 2019_08_21
myY_orignal <- read.table("/Users/zhen0001/Box\ Sync/Pheno_Geno_EBV_SSF_budset3_disease_Northphenology_Sävar_Frost_FFFF.txt", header = TRUE, sep="\t") # 2019_08_21

#********************2019-10-01 *************************
# 1) to select the traits which I used and then do so Pearson correlations 
# 2) do calculate the heritablity  or repeatability, G heritability and 
names(myY_orignal)[1]<-"Taxa"
myY_orignal<-myY_orignal[order(myY_orignal$Taxa),]  # order the sample ids based on Taxa. 
head(myY_orignal); dim(myY_orignal) ; str(myY_orignal)
names(myY_orignal)
table(factor(myY_orignal$Group_FFF))
myY_orignal<-myY_orignal[order(myY_orignal$Taxa),]

trait.myY_orignal<-myY_orignal[!is.na(myY_orignal$Ajd_GV_Pilodyn),]  # order the sample ids based on Taxa. 
trait.myY_orignal<-myY_orignal[!is.na(myY_orignal$Ajd_GV_Velocity),]  # order the sample ids based on Taxa. 
trait.myY_orignal<-myY_orignal[!is.na(myY_orignal$Ajd_GV_MOE),]  # order the sample ids based on Taxa. 
trait.myY_orignal<-myY_orignal[!is.na(myY_orignal$Adj_Skog91_Chen_V2),]  # order the sample ids based on Taxa. 
trait.myY_orignal<-myY_orignal[!is.na(myY_orignal$Adj_Hjd_Scale_DEBV),]  # order the sample ids based on Taxa. 
trait.myY_orignal<-myY_orignal[!is.na(myY_orignal$Adj_DBH_Scale_DEBV),]  # order the sample ids based on Taxa. 
trait.myY_orignal<-myY_orignal[!is.na(myY_orignal$S1150_FD_SDEBV),]  # order the sample ids based on Taxa. 
trait.myY_orignal<-myY_orignal[!is.na(myY_orignal$F1215_FD7_SDEBV),]  # order the sample ids based on Taxa. 
trait.myY_orignal<-myY_orignal[!is.na(myY_orignal$F1146_Frost_SDEBV),]  # order the sample ids based on Taxa. 
trait.myY_orignal<-myY_orignal[!is.na(myY_orignal$F11501215_FD_SDEBV),]  # order the sample ids based on Taxa. 

table(factor(trait.myY_orignal$Group_FFF))


#ALP    C_SE     CEU  Indigo   Kirov     NFE     NPL obovata Omorika     ROM Rus_Bal 
#1125    1571     508      23      16    1041     197      53      21     124     490 
#sum(c(1125  ,  1571  ,   508 ,     23    ,  16 ,   1041   ,  197   ,   53   ,   21  ,   124   ,  490 ))


F1150_adj_Sko91<- drop.levels(subset(myY_orignal, !is.na(myY_orignal$S1150_FD_SDEBV) |!is.na(myY_orignal$F1215_FD7_SDEBV)   ) ) # done  budburst F1146
dim(F1150_adj_Sko91)

# 2020-05-10
ROM_myY_original<-drop.levels(subset(trait.myY_orignal, trait.myY_orignal$Group_FFF=="ROM") ) # done
dim(ROM_myY_original)

ALP_myY_original<-drop.levels(subset(myY_orignal, myY_orignal$Group_FFF=="ALP"))  # done
CSE_myY_original<-drop.levels(subset(myY_orignal, myY_orignal$Group_FFF=="C_SE")) # done
CEU_myY_original<-drop.levels(subset(myY_orignal, myY_orignal$Group_FFF=="CEU") ) # done
NFE_myY_original<-drop.levels(subset(myY_orignal, myY_orignal$Group_FFF=="NFE") ) # done
NPL_myY_original<-drop.levels(subset(myY_orignal, myY_orignal$Group_FFF=="NPL") ) # done
ROM_myY_original<-drop.levels(subset(myY_orignal, myY_orignal$Group_FFF=="ROM") ) # done
Rus_myY_original<-drop.levels(subset(myY_orignal, myY_orignal$Group_FFF=="Rus_Bal") ) # done
Pascal_myY_original<-drop.levels(subset(myY_orignal, !is.na(myY_orignal$Budburst.Observed) ))
S1150_FD_SDEBV<-drop.levels(subset(myY_orignal, !is.na(myY_orignal$S1150_FD_SDEBV)) ) # done 
F1215_FD7_SDEBV<-drop.levels(subset(myY_orignal, !is.na(myY_orignal$F1215_FD7_SDEBV)) ) # done 
budburst_JUN<- drop.levels(subset(myY_orignal, !is.na(myY_orignal$budburst_JUN)) ) # done 
Adj_Budflush17_Savar<- drop.levels(subset(myY_orignal, !is.na(myY_orignal$Adj_Budflush17_Savar)) ) # done 
Adj_Skog91_Ekebo<-drop.levels(subset(myY_orignal, !is.na(myY_orignal$Adj_Skog91_Ekebo)) ) # done 
Adj_wSK091_Ekebo<- drop.levels(subset(myY_orignal, !is.na(myY_orignal$adj_wSK091_Ekebo)) ) # done 
budset_Days2<- drop.levels(subset(myY_orignal, !is.na(myY_orignal$budset_Days2)) ) # done 
Sfrost_7_F1146_SDEBV<- drop.levels(subset(myY_orignal, !is.na(myY_orignal$Sfrost_7_F1146_SDEBV)) ) # done 
Adj_Pilodyn_Savar<- drop.levels(subset(myY_orignal, !is.na(myY_orignal$Adj_Pilodyn_Savar)) ) # done 
Adj_Velocity_Savar<- drop.levels(subset(myY_orignal, !is.na(myY_orignal$Adj_Velocity_Savar)) ) # done 
Adj_MOE_Savar<- drop.levels(subset(myY_orignal, !is.na(myY_orignal$Adj_MOE_Savar)) ) # done 
F1146_Frost_SDEBV<- drop.levels(subset(myY_orignal, !is.na(myY_orignal$F1146_Frost_SDEBV)) ) # done 
F1146_adj_Sko91<- drop.levels(subset(myY_orignal, !is.na(myY_orignal$F1146_Frost_SDEBV)) ) # done  budburst F1146
F1150_adj_Sko91<- drop.levels(subset(myY_orignal, !is.na(myY_orignal$S1150_FD_SDEBV)) ) # done  budburst F1146


Adj_Pilodyn_Brunsberg<- drop.levels(subset(myY_orignal, !is.na(myY_orignal$Adj_Pilodyn_Brunsberg)) ) # done 
Adj_Velocity_Brunsberg<- drop.levels(subset(myY_orignal, !is.na(myY_orignal$Adj_Velocity_Brunsberg)) ) # done 
Adj_MOE_Brunsberg<- drop.levels(subset(myY_orignal, !is.na(myY_orignal$Adj_MOE_Brunsberg)) ) # done 
SWG_tot<- drop.levels(subset(myY_orignal, !is.na(myY_orignal$SWG_tot)) ) # done 
Log_Lesion<- drop.levels(subset(myY_orignal, !is.na(myY_orignal$Log_Lesion)) ) # done 

F11501215_FD<-myY_orignal[!is.na(myY_orignal$F11501215_FD_SDEBV),]


dim(ALP_myY_original) # 1125
dim(CSE_myY_original) # 1571
dim(CEU_myY_original) # 508
dim(NFE_myY_original) # 1241
dim(NPL_myY_original) # 197
dim(ROM_myY_original) # 124
dim(Rus_myY_original) # 490

dim(S1150_FD_SDEBV) # 1151
dim(F1215_FD7_SDEBV)  
dim(budburst_JUN)
dim(Adj_Budflush17_Savar)
dim(Adj_wSK091_Ekebo)
dim(budset_Days2)
dim(Sfrost_7_F1146_SDEBV)
dim(Adj_Skog91_Ekebo)
#dim(Pascal_myY_original) # ca. 800 done
dim(F1146_adj_Sko91)
dim(F1150_adj_Sko91)
dim(F1146_Frost_SDEBV)
dim(F11501215_FD)

# test for specific traits
# myY_orignal_Trait<-droplevels(myY_orignal[!is.na(myY_orignal$Adj_Skog91_Chen),])

# 2020 -5-11
myGD_Trait<-myGD[is.element(myGD$Taxa, ROM_myY_original$Taxa) ,]
dim(myGD_Trait)

#  subset of group
myGD_Trait<-myGD[myY_orignal$Group_FFF=="ALP" & !is.na(myY_orignal$Group_FFF) ,]
myGD_Trait<-myGD[myY_orignal$Group_FFF=="C_SE" & !is.na(myY_orignal$Group_FFF) ,]
myGD_Trait<-myGD[myY_orignal$Group_FFF=="CEU" & !is.na(myY_orignal$Group_FFF) ,]
myGD_Trait<-myGD[myY_orignal$Group_FFF=="NFE" & !is.na(myY_orignal$Group_FFF) ,]
myGD_Trait<-myGD[myY_orignal$Group_FFF=="NPL" & !is.na(myY_orignal$Group_FFF) ,]
myGD_Trait<-myGD[myY_orignal$Group_FFF=="ROM" & !is.na(myY_orignal$Group_FFF) ,]
myGD_Trait<-myGD[myY_orignal$Group_FFF=="Rus_Bal" & !is.na(myY_orignal$Group_FFF) ,]
myGD_Trait<-myGD[!is.na(myY_orignal$Adj_Budflush17_Savar),]
myGD_Trait<-myGD[!is.na(myY_orignal$Adj_Skog91_Ekebo),]

myGD_Trait<-myGD[!is.na(myY_orignal$Adj_Pilodyn_Savar),]
myGD_Trait<-myGD[!is.na(myY_orignal$Adj_Velocity_Savar),]
myGD_Trait<-myGD[!is.na(myY_orignal$Adj_MOE_Savar),]

myGD_Trait<-myGD[!is.na(myY_orignal$Adj_Pilodyn_Brunsberg),]
myGD_Trait<-myGD[!is.na(myY_orignal$Adj_Velocity_Brunsberg),]
myGD_Trait<-myGD[!is.na(myY_orignal$Adj_MOE_Brunsberg),]
myGD_Trait<-myGD[!is.na(myY_orignal$SWG_tot),]
myGD_Trait<-myGD[!is.na(myY_orignal$Log_Lesion),]
myGD_Trait<-myGD[!is.na(myY_orignal$Sfrost_7_F1146_SDEBV),]
myGD_Trait<-myGD[!is.na(myY_orignal$F1146_Frost_SDEBV),]
myGD_Trait<-myGD[!is.na(myY_orignal$F1146_Frost_SDEBV),]

myGD_Trait<-myGD[!is.na(myY_orignal$F1215_FD7_SDEBV),]
myGD_Trait<-myGD[!is.na(myY_orignal$S1150_FD_SDEBV),]

myGD_Trait<-myGD[!is.na(myY_orignal$F11501215_FD_SDEBV),]



# myGD_Trait<-myGD
dim(myGD_Trait)
names(myY_orignal)
myY<-SWG_tot[, c(1, 43)] # disease
myY<-Log_Lesion[, c(1, 44)] # disease
myY<-Sfrost_7_F1146_SDEBV[, c(1, 69)] # Sfrost_7_F1146_SDEBV
myY<-F1146_Frost_SDEBV[, c(1, 122)] # F1146_Frost_SDEBV
myY<-F1146_Frost_SDEBV[, c(1, 123)] # F1146_Hjd_7_SDEBV
myY<-F1146_Frost_SDEBV[, c(1, 124)] # F1146_Dia_12_SDEBV

myY<-F1215_FD7_SDEBV[, c(1, 115)] # F1215_Frost_SDEBV
myY<-F1146_adj_Sko91[, c(1, 9)] # F1146_adj_Sko91
myY<-F1150_adj_Sko91[, c(1, 9)] # F1150_adj_Sko91
myY<-F11501215_FD[, c(1, 140)] #F11501215_FD_SDEBV

head(myY)


# Adj_Skog91_Chen_V2 , bud burst
myY<-ALP_myY_original[, c(1, 73)]
myY<-CSE_myY_original[, c(1, 73)]
myY<-CEU_myY_original[, c(1, 73)] 
myY<-NFE_myY_original[, c(1, 73)] 
myY<-NPL_myY_original[, c(1, 73)]
myY<-ROM_myY_original[, c(1, 73)]
myY<-Rus_myY_original[, c(1, 73)]


#  Adj_Hjd_Scale_DEBV
myY<-ALP_myY_original[, c(1, 23)]
myY<-CSE_myY_original[, c(1, 23)]
myY<-CEU_myY_original[, c(1, 23)] 
myY<-NFE_myY_original[, c(1, 23)] 
myY<-NPL_myY_original[, c(1, 23)]
myY<-ROM_myY_original[, c(1, 23)]
myY<-Rus_myY_original[, c(1, 23)]

#  Adj_DBH_Scale_DEBV
myY<-ALP_myY_original[, c(1, 26)] 
myY<-CSE_myY_original[, c(1, 26)] 
myY<-CEU_myY_original[, c(1, 26)] 
myY<-NFE_myY_original[, c(1, 26)] 
myY<-NPL_myY_original[, c(1, 26)] 
myY<-ROM_myY_original[, c(1, 26)]
myY<-Rus_myY_original[, c(1, 26)]


# Velocity subgroup
myY<-ALP_myY_original[, c(1, 18)] 
myY<-CSE_myY_original[, c(1, 18)]  
myY<-CEU_myY_original[, c(1, 18)] 
myY<-NFE_myY_original[, c(1, 18)] 
myY<-NPL_myY_original[, c(1, 18)] 
myY<-ROM_myY_original[, c(1, 18)]
myY<-Rus_myY_original[, c(1, 18)]

# Pilodyn subgroup
myY<-ALP_myY_original[, c(1, 19)] 
myY<-CSE_myY_original[, c(1, 19)]  
myY<-CEU_myY_original[, c(1, 19)] 
myY<-NFE_myY_original[, c(1, 19)] 
myY<-NPL_myY_original[, c(1, 19)] 
myY<-ROM_myY_original[, c(1, 19)]
myY<-Rus_myY_original[, c(1, 19)]

# MOE subgroup
myY<-ALP_myY_original[, c(1, 20)] 
myY<-CSE_myY_original[, c(1, 20)]  
myY<-CEU_myY_original[, c(1, 20)] 
myY<-NFE_myY_original[, c(1, 20)] 
myY<-NPL_myY_original[, c(1, 20)] 
myY<-ROM_myY_original[, c(1, 20)]
myY<-Rus_myY_original[, c(1, 20)]

# Ekebo

myY<-Adj_Skog91_Ekebo[, c(1, 16)] # budburst

# Savar
myY<-Adj_Budflush17_Savar[, c(1, 27)] # budburst
myY<-Adj_Pilodyn_Savar[, c(1, 28)]
myY<-Adj_Velocity_Savar[, c(1, 29)]
myY<-Adj_MOE_Savar[, c(1, 30)]

# Brunsberg
myY<-Adj_Velocity_Brunsberg[, c(1, 77)] # Velocity
myY<-Adj_Pilodyn_Brunsberg[, c(1, 78)] # Pilodyn
myY<-Adj_MOE_Brunsberg[, c(1, 79)] # start

exclude_Sample<-read.xls("/Users/zhen0001/Box\ Sync/GWAS/results/correct_population_cluster.xlsx", sheet=1)
head(exclude_Sample); dim(exclude_Sample)
sub.exclde_Sample<-drop.levels(subset(exclude_Sample, exclude_Sample$original=="ALP")) # start
sub.exclde_Sample<-drop.levels(subset(exclude_Sample, exclude_Sample$original=="CSE_Height")) # done
sub.exclde_Sample<-drop.levels(subset(exclude_Sample, exclude_Sample$original=="CEU")) # done 
sub.exclde_Sample<-drop.levels(subset(exclude_Sample, exclude_Sample$original=="NPL")) # done 
sub.exclde_Sample<-drop.levels(subset(exclude_Sample, exclude_Sample$original=="ROM")) # start
sub.exclde_Sample<-drop.levels(subset(exclude_Sample, exclude_Sample$original=="RUA_Bal")) # start
sub.exclde_Sample<-drop.levels(subset(exclude_Sample, exclude_Sample$original=="Brunsberg")) # start
sub.exclde_Sample<-drop.levels(subset(exclude_Sample, exclude_Sample$original=="Disease")) # start

myY=subset(myY, !is.element(myY$Taxa, sub.exclde_Sample$taxa ))  
myGD_Trait=subset(myGD_Trait,!is.element(myGD_Trait$taxa, sub.exclde_Sample$taxa )) 

myY[is.na(myY)]=NaN # Gapit need NaN as missing value
head(myY); dim(myY)
dim(myGD_Trait)

#*******only BLINK needed**** this very important to match the same dimension of myGD and  myGM 
myGD_Trait<-t(myGD_Trait[, -1])

dim(myY); dim(myGD_Trait); dim(myGM)
myBlink=Blink(Y=myY,GD=myGD_Trait,GM=myGM, maxLoop = 10, time.cal = T) # Pilodyn, velocity, MoE should not use PC, Pilodyn is not significant if we use PC

#*************************************
# this method is quite slow actually
# ************************************
myGAPIT<- GAPIT(
  Y=myY,
  GD=myGD_Trait,
  GM=myGM,
  model=c("CMLM") # ,   # FarmCPU
  # model=c("BLINK"), 
  # model=c("FarmCPU"),  
  # CV = myCV,
 #  PCA.total=3
  # Model.selection = TRUE,
  # SNP.MAF=0.05,
  # file.output=T
  )
# to get 3D
if(!require(scatterplot3d)) 
  install.packages("scatterplot3d")

