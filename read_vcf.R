library("readxl")
library(TwoSampleMR)
library(phenoscanner)
library(LDlinkR)
library(dplyr)
library(cause)
library(vcfR)
library(stringr)
library(tidyverse)
rm(list=ls())
gc()


server<-'home'
setwd('.../hypertension/data')
vcf_l<-list.files('.../hypertension/data',pattern='*.vcf$',full.names=TRUE,recursive=TRUE)



list_fun<-function(vcf_list){
  
  tmp<-strsplit(basename(vcf_list),split=".",fixed=TRUE)
  b<-unlist(lapply(tmp,head,1))
  vcf_list<-data.frame(vcf_list,b)
  names(vcf_list)<-c('path','id')
  vcf_list$folder<-dirname(vcf_list$path)
  return(vcf_list)
}
vcf_list<-list_fun(vcf_l)

fin_list<-list_fun('.../finngen_R5_R18_ELEVATED_ERYTHROCYTE_SEDIM_RATE_ABNORMALITY_PLASMA_VISCO')
#############################################################################################################################
#preparation for cause


vcf2data<-function(path){
  vcf<-vcfR::read.vcfR(path)
  arGt<-as.data.frame(vcf@gt)
  # arGt<-head(arGt,10)
  # colnames(arGt)<-'GT'
  # arGt[,c(2)]
  colnames(arGt)<-c('colnam1','colnam2')
  # tt1<-head(vcf@fix,10)
  # 
  tt0<-head(vcf@gt,10)
  tt0
  nam<-strsplit(as.character(tt0[1]),split=':')[[1]]
  len<-length(strsplit(as.character(tt0[1]),split=':')[[1]])
  l<-list()
  for (i in 1:length(nam)){
    
    h<-unlist(strsplit(as.character(arGt$colnam2),split=':'))[seq(i,dim(arGt)[1]*len,len)]
    l[[nam[i]]]<-h
    l<-as.data.frame(l)
  }
  l$pval<-10^(-as.numeric(l$LP))
  l<-transform(l, SE = as.numeric(SE),  ES = as.numeric(ES),pval = as.numeric(pval))
  l<-dplyr::rename(l, ID0 = ID)#new_name = old_name
  #l$ID<-NULL
  l<-cbind(as.data.frame(l),as.data.frame(vcf@fix))
  l <- l[!(is.na(l$SE)), ]
  l <- l[!(is.na(l$ES)), ]
  l<-l[l$SE>0,]
  l <- l[!(is.na(l$ID)), ]
  l <- l[l$ID!="", ]
  l<-l[startsWith(l$ID, 'rs'),]#
  l<-l[str_count(l[,c('ID')], c("_rs"))==0,]
  #l$pval<-pnorm(abs(l$ES/l$SE),lower.tail = F)*2
  l<-l%>%subset(pval>0)
  return(l)
}
fin_fun<-function(path){
  fin_cols<-c('chrom','pos','ref','alt','rsids','nearest_genes','pval','beta','sebeta','maf','maf_cases','maf_controls','n_hom_cases','n_het_cases','n_hom_controls','n_het_controls')
  d2<-read.table(path,sep='\t',col.names=fin_cols)
  d2<-as_tibble(d2)
  d2<-d2%>%transform(., pval = as.numeric(pval),  beta = as.numeric(beta),sebeta = as.numeric(sebeta),maf = as.numeric(maf))%>%filter(pval>0,sebeta>0,rsids!="")
  return(d2)
}
ptm <- proc.time()
for (j in 1:nrow(vcf_list)){
  path<-vcf_list$path[j]
  id<-vcf_list$id[j]
  folder<-vcf_list$folder[j]
  # tsv_check<-gsub(".vcf", ".tsv", vcf_list$vcf_list[j], fixed = TRUE)  
  # if (file.exists(tsv_check)){
  #   print(paste0('This ','vcf ',vcf_list$id[j], ' had been done before!'))
  # } else{
  print(paste0('The ',j,' trial!!!'))
  l<-vcf2data(path)
  save(l,file=paste(folder,paste0(id,'.RData'),sep='/'))
  print(head(l))
  rm(l)
  gc()}
#}  
for (j in 1:nrow(fin_list)){
  path<-strsplit(fin_list$path[j],split=".",fixed=TRUE)[[1]][1]
  id<-fin_list$id[j]
  # tsv_check<-gsub(".vcf", ".tsv", vcf_list$vcf_list[j], fixed = TRUE)  
  # if (file.exists(tsv_check)){
  #   print(paste0('This ','vcf ',vcf_list$id[j], ' had been done before!'))
  # } else{
  print(paste0('The ',j,' trial!!!'))
  l<-fin_fun(path)
  save(l,file=paste0(path,'.RData'))
  print(head(l))
  rm(l)
  gc()}
print(paste0('time spent:',' ',proc.time()-ptm))