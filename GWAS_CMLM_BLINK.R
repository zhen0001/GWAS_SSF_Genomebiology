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

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(c("multtest")) # this is only way to install package from BiocManager

library(bigmemory)
library(biganalytics)
library(multtest)
# library(gplots)
library(LDheatmap)
library(genetics)
library(ape)
library(EMMREML)
library(compiler) #this library is already installed in R
library("scatterplot3d")
source("http://zzlab.net/GAPIT/gapit_functions.txt")


# The EMMA library was developed by Kang et al. (2008). One line was added to the library to handle the
# grid search with "NA" likelihood. The modified library can be installed by typing this command line:
source("http://zzlab.net/GAPIT/emma.txt")
setwd("/Users/zhen0001/Documents/onedrive/BLINK-5238/")

require(data.table) # fread is similar to read.table, but faster and convenient all controls suchas sep, colcloasses, nrows,
data<-fread("out.012", na.strings=c("",".", "NA")) # "-1" it is great helpful for # mdata1[mdata1=="-1"]=NA
data<-data[,-1]
indv<-fread("out.012.indv", header = FALSE) # the name of individuals
sortmarker_id<- fread("out.012.pos", header = FALSE) # the name of markers
length(unique(sortmarker_id$V1)) # to count the number of contig

colnames(data)<-paste(sortmarker_id$V1, sortmarker_id$V2, sep="_") # [1:10000]
row.names(data)<-indv$V1
myGD<-data.frame(indv$V1, data)
colnames(myGD)[1]="Taxa"
SNP<-paste(sortmarker_id$V1, sortmarker_id$V2, sep="_")
myGM<-data.frame(SNP,c(rep(1:10, each=29000),rep(10, 9240)) , 1:length(SNP))
colnames(myGM)<-c("rs", "chr", "pos") # rs: reference snp
rm(SNP); rm(data)

#*****************************************************
# BLINK in R version
# association analysis
source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("http://zzlab.net/FarmCPU/FarmCPU_functions.txt")
library(BLINK)


# to get the frequency of each SNPs
freq<-function(x){
  (sum(x==2)*2+ sum(x==1))/(length(x)*2)
}

freq.snp<-apply(myGD[, -1],2, freq)
myGD<-myGD[, c(TRUE, freq.snp>=0.03 & freq.snp<=0.97)]
myGM<-myGM[c(freq.snp>=0.03 & freq.snp<=0.97), ]

# GWAS phenotype data
myY_orignal <- read.table("/Users/zhen0001/Box\ Sync/GB_backup/Zendo_databacp/Phenotypic_data_GWAS_5056.txt", header = TRUE, sep="\t") # 2019_08_21
names(myY_orignal)[1]<-"Taxa"
myY_orignal<-myY_orignal[order(myY_orignal$Taxa),]  # order the sample ids based on Taxa. 
table(factor(myY_orignal$Group_FFF))

myGD<-myGD[is.element(myGD$Taxa, myY_orignal$Taxa),]

# Matrix transpose 
T.myGD<-t(myGD[, -1])

# BLINK without PC
for (i in 3:9) {
  myY<-myY_orignal[,c(1,i )]
  myBlink=Blink(Y=myY, GD=T.myGD, GM=myGM, maxLoop = 10, time.cal = T) 
}

#**************************
# a PC to a few PCs
#**************************
library(AGHmatrix)
G<-Gmatrix(as.matrix(myGD[,-1]), ploidy=2, missingValue=-1, maf=0.03, method="VanRaden")
pca = prcomp(as.matrix(G), center = TRUE,scale. = TRUE)
pc1<-as.matrix(pca$x[,1])

for (i in 3:9) {
  myY<-myY_orignal[, c(1,i)]
  myBlink=Blink(Y=myY, GD=myGD, GM=myGM, CV = pc1, maxLoop = 10, time.cal = T) 
}

#************************
# Gapit CMLM
#**************************
# wthout PC
for (i in 3:9) {
myY<-myY_orignal[, c(1,i)]
  
myGAPIT<- GAPIT(
  Y=myY,
  GD=myGD,
  GM=myGM,
  model=c("CMLM") )
}

# with PC 
for (i in 3:9) {
  myY<-myY_orignal[, c(1,i)]
  myGAPIT<- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    model=c("CMLM"),
    PCA.total=3)
}



