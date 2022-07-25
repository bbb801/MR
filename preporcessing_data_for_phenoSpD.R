
library("readxl")
library(TwoSampleMR)
library(phenoscanner)
library(LDlinkR)
library(dplyr)
library(cause)
library(RadialMR)
library(MVMR)
library(MRPRESSO)
library(boot.pval)
library(tidyverse)
rm(list=ls())
gc()



server<-'home'
data_list<-list.files('.../2_20/data/figgen',pattern='.RData$',full.names=TRUE,recursive=TRUE)
whole_result_folder<-'.../2_20/data/summary_data_for_corr'
data_table<-'.../2_20/data_table.csv'
source('.../2_20/code/umr_tool.R')

dir.create(whole_result_folder, showWarnings = TRUE, recursive = TRUE)
setwd(whole_result_folder)



id_exposure <- c("ukb-d-30040_irnt")
id_outcome<-c( "finn-b-ST19_TRAUMAT_SUBAR_HAEMORRHAGE", 
               "finn-b-I9_SAH",  "finn-b-I9_INTRACRA", 
               "finn-b-I9_ICH")
MCV_vcf<-data_list[str_detect(data_list,id_exposure)]


final<-d_fun(data_list,data_table)%>%arrange(trait)
final<-final[c(4,6,2,7,8,5,1,3),]
final<-final%>%filter(.,!trait %in% c("Traumatic subarachnoid haemorrhage","Intracranial trauma"))

dput(final$trait)
dput(final$data_list)
final$trait2<-c("MCV","RCDW","ICH", "SAH",
  "nITH",  
  "AOS")

ptm <- proc.time()
for (i in 1:nrow(final)){
  if (i==1){
    load(final$data_list[i])
    l<-l[startsWith(l$ID, 'rs'),]#
    l<-l[str_count(l[,c('ID')], c("_rs"))==0,]
    l1<-l%>%as_tibble%>%rename(MCV_b=ES,MCV_se=SE,allele_0=REF,allele_1=ALT)%>% distinct(ID, .keep_all = TRUE)%>%drop_na(ID)%>%filter(.,0<AF,pval>0,AF<1)%>%select(allele_0,allele_1,MCV_b,MCV_se,ID)
    rm(l)
    gc()
  } else if (i==2){
    load(final$data_list[i])
    l<-l[startsWith(l$ID, 'rs'),]#
    l<-l[str_count(l[,c('ID')], c("_rs"))==0,]
    l<-l%>%as_tibble%>%rename(RCDW_b=ES,RCDW_se=SE)%>% distinct(ID, .keep_all = TRUE)%>%drop_na(ID)%>%filter(.,0<AF,pval>0,AF<1)%>%select(RCDW_b,RCDW_se,ID)%>%inner_join(l1, 
                                                                                                                                                                          c("ID" = "ID"))
    rm(l1)
    l2<-l
    rm(l)
    gc()
  } else if (i==3) {
    load(final$data_list[i])
    l<-l[startsWith(l$rsids, 'rs'),]#
    l<-l[str_count(l[,c('rsids')], c("_rs"))==0,]
    l<-l%>%as_tibble%>%rename(ICH_b=beta,ICH_se=sebeta,ID=rsids)%>% distinct(ID, .keep_all = TRUE)%>%drop_na(ID)%>%filter(.,0<maf,pval>0,maf<1)%>%select(ICH_b,ICH_se,ID)%>%inner_join(l2, 
                                                                                                                                                                                   c("ID" = "ID"))
    rm(l2)
    l3<-l
    rm(l)
    gc()
  } else if (i==4) {
    load(final$data_list[i])
    l<-l[startsWith(l$rsids, 'rs'),]#
    l<-l[str_count(l[,c('rsids')], c("_rs"))==0,]
    l<-l%>%as_tibble%>%rename(SAH_b=beta,SAH_se=sebeta,ID=rsids)%>% distinct(ID, .keep_all = TRUE)%>%drop_na(ID)%>%filter(.,0<maf,pval>0,maf<1)%>%select(SAH_b,SAH_se,ID)%>%inner_join(l3, 
                                                                                                                                                                                   c("ID" = "ID"))
    rm(l3)
    l4<-l
    rm(l)
    gc()
  } else if (i==5) {
    load(final$data_list[i])
    l<-l[startsWith(l$rsids, 'rs'),]#
    l<-l[str_count(l[,c('rsids')], c("_rs"))==0,]
    l<-l%>%as_tibble%>%rename(nITH_b=beta,nITH_se=sebeta,ID=rsids)%>% distinct(ID, .keep_all = TRUE)%>%drop_na(ID)%>%filter(.,0<maf,pval>0,maf<1)%>%select(nITH_b,nITH_se,ID)%>%inner_join(l4, 
                                                                                                                                                                                       c("ID" = "ID"))
    rm(l4)
    l5<-l
    rm(l)
    gc()
  } else if (i==6) {
    load(final$data_list[i])
    l<-l[startsWith(l$rsids, 'rs'),]#
    l<-l[str_count(l[,c('rsids')], c("_rs"))==0,]
    l<-l%>%as_tibble%>%rename(AOS_b=beta,AOS_se=sebeta,ID=rsids)%>% distinct(ID, .keep_all = TRUE)%>%drop_na(ID)%>%filter(.,0<maf,pval>0,maf<1)%>%select(AOS_b,AOS_se,ID)%>%inner_join(l5, 
                                                                                                                                                                                   c("ID" = "ID"))
    rm(l5)
    l<-as.data.frame(l)
    rownames(l)<-(l$ID)
   # l<-l%>%select(-ID)
    
    gc()
  }
}
dput(names(l))
l<-l[,c("allele_0", "allele_1", "MCV_b", "MCV_se", "RCDW_b", 
  "RCDW_se","ICH_b", "ICH_se", "SAH_b", "SAH_se", "nITH_b", "nITH_se", "AOS_b", "AOS_se"
   )]
write.table(l, file =paste(whole_result_folder,'_4_traits_data_for_corr.txt', sep ="/"), row.names =TRUE, col.names =TRUE, quote =FALSE)
print(paste0('time spent:',' ',proc.time()-ptm))