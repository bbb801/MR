library("readr")
library("tidyverse")
library(coloc)
library(tidyverse)
library(dplyr)
library(stringr)
library(GWAS.utils)
rm(list=ls())
gc()
clump_r2 = 0.01
maf_threshold = 0.01
presso_needed=FALSE

server<-'home'
har_list<-list.files('.../2_20/result_save/umr_clump_0.01/ukb-d-30040_irnt',pattern='final_har.RData$',full.names=TRUE,recursive=TRUE)
data_list<-list.files('.../2_20/data/figgen',pattern='*.RData$',full.names=TRUE,recursive=TRUE)
whole_folder<-paste0('.../2_20/result/coloc_clump_',clump_r2)
data_table<-'.../2_20/data_table.csv'
source('.../2_20/code/umr_tool.R')
bfile<-".../plink_reference/EUR"

print(rep(1,20))

#mr_package_input_snp_corr<-ld_matrix(raddat$SNP, with_alleles = TRUE, pop = "EUR")
#unlink(whole_folder,recursive = TRUE)
dir.create(whole_folder, showWarnings = TRUE, recursive = TRUE)
ptm <- proc.time()

#source(mvmr_tool)

id_exposure <- c("ukb-d-30040_irnt")
id_outcome<-c( "finn-b-ST19_TRAUMAT_SUBAR_HAEMORRHAGE", 
               "finn-b-I9_SAH",  "finn-b-I9_INTRACRA", 
               "finn-b-I9_ICH")
MCV_vcf<-data_list[str_detect(data_list,id_exposure)]


final<-d_fun(data_list,data_table,id=c(id_outcome))

har_list<-data.frame(har_list)
har_list$id<-unlist(lapply(strsplit(har_list$har_list,split='/'),path_fun1))

final<-inner_join(final, har_list, 
                  c("id" = "id"))%>%rename(out_vcf=data_list)
final$exp_vcf<-MCV_vcf
rm(har_list)
gc()
#i<-1

##################################

RCDW_vcf<-data_list[str_detect(data_list,'ebi-a-GCST006804')]
RCDW_final<-final[final$id=="finn-b-ST19_TRAUMAT_SUBAR_HAEMORRHAGE",]
RCDW_final$exp_vcf<-RCDW_vcf
RCDW_final$har_list<-NA
final<-rbind(final,RCDW_final)
##########################

ptm <- proc.time()
exp_id<-"ukb-d-30040_irnt"
#exp_id<-'ebi-a-GCST006804'###############################################################################
if (exp_id=="ukb-d-30040_irnt"){
  exp_N<-350473
  exp_trait<-'MCV'
  loop<-1:length(final$out_vcf)
} else if (exp_id=='ebi-a-GCST006804'){
  exp_N<-116666
  exp_trait<-'RCDW'
  loop<-length(final$out_vcf):length(final$out_vcf)
}
#for (i in loop){print(i)}



for (i in loop){#1:length(final$data_list)
  
  print(paste('The ', i, ' trial！'))
  
  #  re_log_folder<-'radial_outliner.log'
  #  unlink(re_log_folder,recursive = TRUE)
  #  dir.create(f, showWarnings = TRUE, recursive = TRUE)
  #  con <- file(re_log_folder) # 
  #  sink(con, append=TRUE) # 
  #  sink(con, append=TRUE, type="message") # 
  load(final$out_vcf[i])
  
  vcf_dat_out<-l%>% distinct(rsids, .keep_all = TRUE)%>%filter(.,maf<1,0<maf,pval>0)%>%drop_na(rsids)
  vcf_dat_out<-vcf_dat_out[startsWith(vcf_dat_out$rsids, 'rs'),]#
  vcf_dat_out<-vcf_dat_out[str_count(vcf_dat_out[,c('rsids')], c("_rs"))==0,]
  rm(l)
  load(final$exp_vcf[i])
  vcf_dat_exp<-l%>% distinct(ID, .keep_all = TRUE)%>%transform(., AF=as.numeric(AF))%>%filter(.,0<AF,pval>0,AF<1)%>%transform(., MAF = pmin(AF, 1-AF),POS=as.numeric(POS))
  vcf_dat_exp
  rm(l)
  gc()

  out_id<-final$id[i]
  folder<-paste(whole_folder,exp_id,out_id,sep='/')
  #  unlink(folder,recursive = TRUE)
  dir.create(folder, showWarnings = TRUE, recursive = TRUE)
  setwd(folder)
  
  out_trait<-final$trait[i]
  print(paste('Exp: ',exp_trait,'; ID ',exp_id))
  print(paste('Outcome: ',out_trait,'; ID ',out_id))
  
  ##################
  
  r2_thresh = 0.01
  pval_thresh = 1e-3
  X <- gwas_merge(vcf_dat_exp,vcf_dat_out, snp_name_cols = c("ID","rsids"), 
                  beta_hat_cols = c("ES","beta"), 
                  se_cols = c("SE","sebeta"), 
                  A1_cols = c( "REF","ref"), 
                  A2_cols = c( "ALT","alt"), 
                  pval_cols = c("pval", "pval"))%>%inner_join(vcf_dat_out[,c('rsids','maf','chrom','pos')], 
                                                              c("snp" = "rsids"))%>%inner_join(vcf_dat_exp[,c('ID','MAF')], 
                                                                                               c("snp" = "ID"))%>%transform(.,varbe1=seb1*seb1,varbe2=seb2*seb2)%>%rename(MAF2=MAF,MAF1=maf)
  save(X,file='merge_exp_out.RData')
  rm(vcf_dat_out,vcf_dat_exp)
  gc()
#   clump_path<-paste(folder,paste0('X_clump_r_',r2_thresh,'_p_',pval_thresh,'.RData'),sep='/')
# #  if(!file.exists(clump_path)){
#   X_clump <- X %>%
#     rename(rsid = snp,
#            pval = p1) %>%
#     ieugwasr::ld_clump(dat = .,
#                        clump_r2 = r2_thresh,
#                        clump_p = pval_thresh,
#                        #clump_kb = kb_thresh,
#                        plink_bin = genetics.binaRies::get_plink_binary(), 
#                        bfile = bfile,
#                        pop = "EUR")
#   save(X_clump,file=clump_path) 
# #  } else {
# #    load(clump_path)
# #  }
#   #  rm(X)
#   gc()
#   
#   save(X_clump,file='merge_clump_exp_out.RData')
#   if(server=='hpc'){
#     presso_needed<-TRUE
#   }
#   out_exp_har<-filter_fun(X_clump,folder,out_id=out_id,exp_id=exp_id,out_trait=out_trait,exp_trait=exp_trait,bfile,presso_needed=presso_needed,pattern='merge_clump')
#   save(out_exp_har,file='merge_clump_filter_exp_out.RData')
#   rm(out_exp_har)

  X<-X%>%drop_na(p1,p2)%>%arrange(p1,p2)
  row.names(X)<-NULL
  result_whole <- coloc.abf(dataset1=list(pvalues=X$p1, type="quant", N=exp_N,snp=X$snp), dataset2=list(pvalues=X$p2, type="cc", s=final$ncase[i]/final$sample_size[i],N=final$sample_size[i],snp=X$snp), MAF=X$MAF1)
  save(result_whole,file=paste0('merge_exp_out_whole_MAFexp_abf',nrow(X),'.RData'))
  X_1M<-X%>%filter(.,pos<(X$pos[1]+1000000))%>%filter(.,pos>=(X$pos[1]-1000000))
  row.names(X_1M)<-NULL
  save(X_1M,file='merge_exp_out_1M.RData')  
  X_3M<-X%>%filter(.,pos<(X$pos[1]+3000000))%>%filter(.,pos>=(X$pos[1]-3000000))
  row.names(X_3M)<-NULL
  save(X_3M,file='merge_exp_out_3M.RData')
  X_200K<-X%>%filter(.,pos<(X$pos[1]+200000))%>%filter(.,pos>=(X$pos[1]-200000))
  row.names(X_200K)<-NULL
  save(X_200K,file='merge_exp_out_200K.RData')
  rm(X)
  gc()
  abf_fun<-function(X,d=NA){
    X<-X%>%drop_na(p1,p2)%>%arrange(p1,p2)
    
    if (!is.na(d)){
      print(paste0('Cut regions'))
      X<-X%>%filter(.,pos<(X$pos[1]+d))%>%filter(.,pos>=(X$pos[1]-d))}
    row.names(X)<-NULL
     # save(X_cut,file=paste0('merge_exp_out_',d,'.RData'))  
    result <- coloc.abf(dataset1=list(pvalues=X$p1, type="quant", N=exp_N,snp=X$snp), dataset2=list(pvalues=X$p2, type="cc", s=final$ncase[i]/final$sample_size[i],N=final$sample_size[i],snp=X$snp), MAF=X$MAF1)
      #save(result2,file=paste0('merge_exp_out_d',d,'_MAFexp_abf_snp',nrow(X_cut),'.RData'))
    return(list('X'=X,'result'=result))
  }
  
  abf_fun(X,1000000)
  
  result1 <- coloc.abf(dataset1=list(pvalues=X_1M$p1, type="quant", N=exp_N,snp=X_1M$snp), dataset2=list(pvalues=X_1M$p2, type="cc", s=final$ncase[i]/final$sample_size[i],N=final$sample_size[i],snp=X_1M$snp), MAF=X_1M$MAF1)
  save(result1,file=paste0('merge_exp_out_1M_MAFexp_abf',nrow(X_1M),'.RData'))
  print(paste0('1M abf:',result1))
  rm(result1)
  result2 <- coloc.abf(dataset1=list(pvalues=X_3M$p1, type="quant", N=exp_N,snp=X_3M$snp), dataset2=list(pvalues=X_3M$p2, type="cc", s=final$ncase[i]/final$sample_size[i],N=final$sample_size[i],snp=X_3M$snp), MAF=X_3M$MAF1)
  save(result2,file=paste0('merge_exp_out_3M_MAFexp_abf',nrow(X_3M),'.RData'))
  print(paste0('3M abf:',result2))
  rm(result2)
  result3 <- coloc.abf(dataset1=list(pvalues=X_200K$p1, type="quant", N=exp_N,snp=X_200K$snp), dataset2=list(pvalues=X_200K$p2, type="cc", s=final$ncase[i]/final$sample_size[i],N=final$sample_size[i],snp=X_200K$snp), MAF=X_200K$MAF1)
  save(result2,file=paste0('merge_exp_out_3M_MAFexp_abf',nrow(X_3M),'.RData'))
  print(paste0('3M abf:',result2))
  rm(result2)

  rm(X_1M,X_3M)
  gc()
  if(server=='hpc'){
  try({ld<-ieugwasr::ld_matrix(
    X_1M$snp,
   # X_200K$snp,
    plink_bin = genetics.binaRies::get_plink_binary(),
    bfile = bfile
  )
  save(ld,file='merge_exp_out_1M_ld.RData')})
  rm(ld)
  }

  gc()
}

dataset1<-list(pvaluefinals=X_1M$p1, type="quant", N=exp_N,snp=X_1M$snp)
dataset2<-list(pvalues=X_1M$p2, type="cc", s=final$ncase[i]/final$sample_size[i],N=final$sample_size[i],snp=X_1M$snp,MAF=X_1M$MAF1)
X<-list('dataset1'=dataset1,'dataset2'=dataset2)
library(gassocplot)
library(gwasglue)
temp<-coloc_to_gassocplot(X, bfile = bfile, plink_bin = genetics.binaRies::get_plink_binary())
temp <- coloc_to_gassocplot(result1)
#> Extracting LD matrix for 148 variants
#> Please look at vignettes for options on running this locally if you need to run many instances of this command.
#> Warning in ieugwasr::ld_matrix(markers[["marker"]], with_alleles = FALSE, : The following variants are not present in the LD reference panel
#> rs8110695
#> rs2738459
#> rs17001200
#> Found 145 variants in LD reference panel
gassocplot::stack_assoc_plot(temp$markers, temp$z, temp$corr, traits=temp$traits)
print(paste0('time spent:',' ',proc.time()-ptm))