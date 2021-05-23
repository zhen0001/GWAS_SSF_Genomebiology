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
source("/Users/zhen0001/Box\ Sync/FunctionR/Chen.function.R")

# The EMMA library was developed by Kang et al. (2008). One line was added to the library to handle the
# grid search with "NA" likelihood. The modified library can be installed by typing this command line:

source("http://zzlab.net/GAPIT/emma.txt")
setwd("/Users/zhen0001/Documents/onedrive/BLINK-5238/")

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
myGD[1:10, 1:5]
head(sortmarker_id)
n<-dim(sortmarker_id)[1]

sortmarker_id[sortmarker_id$V1=="MA_12842", ]


SNP<-paste(sortmarker_id$V1, sortmarker_id$V2, sep="_")
myGM<-data.frame(SNP,c(rep(1:10, each=29000),rep(10, 9240)) , 1:n)
colnames(myGM)<-c("rs", "chr", "pos") # rs: reference snp

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
head(freq.snp[order(freq.snp)])

str(freq.snp>=0.03)
sum(freq.snp>=0.03) # 135902
sum(freq.snp>=0.03 & freq.snp<=0.97 ) # 134605 # some snps is not the 2 is minor allele frequency
sum(freq.snp>=0.05 & freq.snp<=0.95 ) # 98187 # so there is something different

# to check if the MyGM$ rs == anno.snp$rs 
anno.snp<-read.table("/Users/zhen0001/Documents/onedrive/BLINK-5238/SSF_indel_Bial_mis0.3.d140_GQ10_mis0.7.maf0.005_HWELocus.Imputed.sorted.P5238.ann.txt", header = TRUE, sep="\t", stringsAsFactors = T)
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

#********************2019-10-01 *************************
# 1) to select the traits which I used and then do so Pearson correlations 
# 2) do calculate the heritablity  or repeatability, G heritability and 
# the whole SSF phenotype data
myY_orignal <- read.table("/Users/zhen0001/Box\ Sync/Pheno_Geno_EBV_SSF_budset3_disease_Northphenology_Sävar_Frost_FFFF.txt", header = TRUE, sep="\t") # 2019_08_21
names(myY_orignal)[1]<-"Taxa"
myY_orignal<-myY_orignal[order(myY_orignal$Taxa),]  # order the sample ids based on Taxa. 
head(myY_orignal); dim(myY_orignal) ; str(myY_orignal)
names(myY_orignal)
table(factor(myY_orignal$Group_FFF))
myY_orignal<-myY_orignal[order(myY_orignal$Taxa),]


myY_orignal<-myY_orignal[, c(1, 4, 6, 73, 23, 26, 19,18, 20)]
head(myY_orignal); dim(myY_orignal)
myY_orignal<-myY_orignal[!is.na(myY_orignal[,4]) | !is.na(myY_orignal$Adj_Hjd_Scale_DEBV) | !is.na(myY_orignal$Adj_DBH_Scale_DEBV) | !is.na(myY_orignal$Ajd_GV_Pilodyn) | !is.na(myY_orignal$Ajd_GV_Velocity) | !is.na(myY_orignal$Ajd_GV_MOE) ,  ]
head(myY_orignal); dim(myY_orignal)

# myY_orignal<-drop.levels(subset(myY_orignal, myY_orignal$Group_FFF=="ALP"))  # done
# myY_orignal<-drop.levels(subset(myY_orignal, myY_orignal$Group_FFF=="C_SE")) # done

# myY_orignal<-drop.levels(subset(myY_orignal, myY_orignal$Group_FFF=="CEU") ) # done
# myY_orignal<-drop.levels(subset(myY_orignal, myY_orignal$Group_FFF=="NFE") ) # done

# myY_orignal<-drop.levels(subset(myY_orignal, myY_orignal$Group_FFF=="NPL") ) # done
#  myY_orignal<-drop.levels(subset(myY_orignal, myY_orignal$Group_FFF=="ROM") ) # done
#  myY_orignal<-drop.levels(subset(myY_orignal, myY_orignal$Group_FFF=="Rus_Bal") ) # done

#*******only BLINK needed**** this very important to match the same dimension of myGD and  myGM 
load("/Users/zhen0001/Documents/onedrive/BLINK-5238/CumulativePVE/G.Rdata")

trait.myGD<-myGD[is.element(myGD$Taxa,myY_orignal$Taxa), ] # here I used 
dim(trait.myGD)
trait.myGD[1:5,1:5]
G[1:5, 1:5]; dim(G)
G<-G[is.element(myGD$Taxa,myY_orignal$Taxa), is.element(myGD$Taxa, myY_orignal$Taxa)]
pca = prcomp(as.matrix(G), center = TRUE,scale. = TRUE)

pc1<-as.matrix(pca$x[,1])
#pc1<-data.frame(Taxa=trait.myGD$Taxa,  PC1=pca$x[,1]) # this is for Gapit
dim(pc1); head(pc1) 
rm(G)

#************************
trait.myGD<-t(trait.myGD[, -1])
dim(myY_orignal); dim(trait.myGD)

for (i in 4:5) {
  myY<-myY_orignal[,c(1,i )]
  myBlink=Blink(Y=myY, GD=trait.myGD, GM=myGM,  CV = pc1, maxLoop = 10, time.cal = T) 
}


dim(myY_orignal); dim(trait.myGD); dim(myGM)
myBlink=Blink(Y=myY_orignal,GD=trait.myGD,GM=myGM, maxLoop = 10, time.cal = T) # Pilodyn, velocity, MoE should not use PC, Pilodyn is not significant if we use PC

#*************************************
# this method is quite slow actually
# ************************************

for (i in 4:9) {
  myY<-myY_orignal[,c(1,i )]
  myGAPIT<- GAPIT(
    Y=myY,
    GD=trait.myGD,
    GM=myGM,
    model=c("CMLM") 
  )}

# optimization for number of PCs based on BIC
myGAPIT<- GAPIT(
  Y=myY,
  GD=trait.myGD,
  GM=myGM,
   model=c("GLM"), 
  #  model=c("CMLM") 
  #  model=c("BLINK"), 
  #   model=c("FarmCPU"),  
  #   CV = pc1
 PCA.total=3
#  Model.selection = TRUE
)
