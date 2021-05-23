
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
library(gplots)
library(LDheatmap)
library(genetics)
library(ape)
library(EMMREML)
library(compiler) #this library is already installed in R
library("scatterplot3d")
library(dplyr)
library(data.table) # fread
source("http://zzlab.net/GAPIT/gapit_functions.txt")

# The EMMA library was developed by Kang et al. (2008). One line was added to the library to handle the
# grid search with "NA" likelihood. The modified library can be installed by typing this command line:
source("http://zzlab.net/GAPIT/emma.txt")

# phenotypc data
mdat <- read.table("/Users/zhen0001/Box\ Sync/GB_backup/Zendo_databacp/S1389_mydata_budburststage.csv" , header = TRUE, sep=",", na.strings = c(".", "", "NA") )
head(mdat)

##################################
# get genotype
#*****************************
# here, we will have 914 genotypes

myGD<-fread("/Users/zhen0001/Box\ Sync/GB_backup/Zendo_databacp/S1389.SNPArrayValidation.txt", na.strings=c("",".", "NA"))
head(myGD[, 1:5])
names(myGD)[1] <- "Taxa"  


Allele.freq<-function(x){
  (sum(x==2)*2+ sum(x==1))/(length(x)*2)
}


# filter some loci with MAF <0.05, in GWAS
freq.snp<-apply(myGD [,-1], 2, Allele.freq) # this is just calcualte allele frequency

# hist(freq.snp); length(freq.snp)
SNPMAF0.05<-names(freq.snp)[freq.snp>=0.05 & freq.snp<=0.95]
myGD<-as.data.frame(myGD)
myGD<-myGD[ , c("Taxa", SNPMAF0.05)] 
myGM<-data.frame(rs=colnames(myGD[,-1]), chr=1, pos=1:dim(myGD[, -1])[2] )


myY_orignal<-  mdat[match( as.character(myGD$Taxa), as.character(mdat$Genotype_id)) ,]
identical(as.character(myY_orignal$Genotype_id),as.character(myGD$Taxa))

# phenotypic data
names(myY_orignal)


  myY<- myY_orignal[, c(2, 1)]
  myY[is.na(myY)]=NaN
  myGAPIT<- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    model=c("Blink"),
    cutOff=0.01 

  )

