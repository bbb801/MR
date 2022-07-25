library(TwoSampleMR)
rm(list=ls())
gc()
library("readxl")
library(TwoSampleMR)
library(tidyverse)
library(RadialMR)
library(GWAS.utils)
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
library(MendelianRandomization)
library(RMVMR)

server<-'home'
corfile<-'.../2_20/result_save/rmvmr_clump_0.01/_corr.txt'
whole_result_folder<-'.../hypertension/systolic'#b1
whole_result_folder<-'.../hypertension/exp_systolic2ich'#b2
cvs_path<-'.../hypertension/alldata.csv'
mvmr_tool<-'.../2_20/code/mvmr_tool.R'
bfile<-"F:/home/featurize/work/gwas_pairs/plink_reference/EUR"
source('.../2_20/code/umr_tool.R')


###setting
source(mvmr_tool)
b<-'b1'
ld_matrix_cal<-FALSE
gencov_cal<-FALSE
develop<-TRUE
train_type<-'mvmr'
ptm <- proc.time()
style<-'Amrican'

d<-read_csv(cvs_path)
sample_rename_fun<-function(d_exp){
  d_exp[is.na(d_exp$sample_size),]$sample_size<-d_exp[is.na(d_exp$sample_size),]$ncase+d_exp[is.na(d_exp$sample_size),]$ncontrol
  d_exp$abbr<-NA
  d_exp$abbr[d_exp$id=='ebi-a-GCST006804']<-'RCDW'
  d_exp$abbr[d_exp$id=='ukb-d-30040_irnt']<-'MCV'
  d_exp$abbr[d_exp$id=='ieu-b-38']<-'SBP'
  d_exp$abbr[d_exp$id=='ieu-b-39']<-'DBP'
  d_exp$abbr[d_exp$id=='finn-b-R18_ELEVATED_ERYTHROCYTE_SEDIM_RATE_ABNORMALITY_PLASMA_VISCO']<-'BD'
  d_exp$abbr[d_exp$id=='finn-b-I9_ICH']<-'ICH'
  d_exp$abbr[d_exp$id=='finn-b-I9_INTRACRA']<-'nITH'
  return(d_exp)
}

out_name<-function(x){
  if(x=='ebi-a-GCST006804'){
    return('RCDW')
  } else if (x=='ukb-d-30040_irnt'){
    return('MCV')
  } else if
    (x=='ieu-b-38'){
      return('SBP')
    }else if
      (x=='ieu-b-39'){
        return('DBP')
        }else if
  (x=='finn-b-R18_ELEVATED_ERYTHROCYTE_SEDIM_RATE_ABNORMALITY_PLASMA_VISCO'){
    return('BD')
      
  }else if
  (x=='finn-b-I9_ICH'){
    return('ICH')
    
  }else if(x=='finn-b-I9_INTRACRA'){
    return('nITH')
  } else{
    return(x)
  }
}

change_name_mvmr<-function(nam,id_outcome,style='American'){
  n<-length(nam)/2
  for (i in 1:length(nam)){
    
    if (i<=n){
      if (nam[i]==out_name(nam[i])){
        outcome_abbr<-out_name_new_list(id_outcome,id2abbr_real=1,reverse=FALSE,style=style)
        nam[i]<-paste0(outcome_abbr,'_b')
      } else{
        nam[i]<-paste0(out_name(nam[i]),'_b')
      }
    } else {
      if (nam[i]==out_name(nam[i])){
        outcome_abbr<-out_name_new_list(id_outcome,id2abbr_real=1,reverse=FALSE,style=style)
        nam[i]<-paste0(outcome_abbr,'_se')
      }else{
        nam[i]<-paste0(out_name(nam[i]),'_se')
      }
    }
  }
  re<-list('renam'=nam,'out_abbr'=outcome_abbr)
  return(re)
}

id_fun<-function(b){
  d<-read_csv(cvs_path)
  exp_id<-c('ebi-a-GCST006804','ukb-d-30040_irnt')
  exp_list<-d%>%filter(id %in%exp_id)
  exp_list$abbr<-out_name_new_list(exp_list$id,id2abbr_real=1,reverse=FALSE,style=style)
  if (b=='b0'){
    out_id<-c("finn-b-I9_SAHANEUR", "finn-b-ST19_TRAUMAT_SUBAR_HAEMORRHAGE", 
    "finn-b-I9_SAH", "finn-b-ST19_INTRACRATRAUMA", "finn-b-I9_INTRACRA", 
    "finn-b-I9_ICH")
    out_list<-d%>%filter(id %in%out_id)
    out_list$abbr<-out_name_new_list(out_list$id,id2abbr_real=1,reverse=FALSE,style=style)
    dat<-list('exp_list'=exp_list,'out_list'=out_list,'exp_id'=exp_id,'out_id'=out_id)
    return(dat)
    
  }else if(b=='b2'){
    med_id<-c("ieu-b-39","ieu-b-38","finn-b-R18_ELEVATED_ERYTHROCYTE_SEDIM_RATE_ABNORMALITY_PLASMA_VISCO")
    med_list<-d%>%filter(id %in%med_id)
    med_list$abbr<-out_name_new_list(med_list$id,id2abbr_real=1,reverse=FALSE,style=style)
    out_id<-c('finn-b-I9_INTRACRA','finn-b-I9_ICH')
    out_list<-d%>%filter(id %in%out_id)
    out_list$abbr<-out_name_new_list(out_list$id,id2abbr_real=1,reverse=FALSE,style=style)
    dat<-list('exp_list'=exp_list,'med_list'=med_list,'out_list'=out_list,'exp_id'=exp_id,'med_id'=med_id,'out_id'=out_id)
  }else if(b=='b1'){
    out_id<-c("ieu-b-39","ieu-b-38","finn-b-R18_ELEVATED_ERYTHROCYTE_SEDIM_RATE_ABNORMALITY_PLASMA_VISCO")
    out_list<-d%>%filter(id %in%out_id)
    out_list$abbr<-out_name_new_list(out_list$id,id2abbr_real=1,reverse=FALSE,style=style)
    dat<-list('exp_list'=exp_list,'out_list'=out_list,'exp_id'=exp_id,'out_id'=out_id)
    return(dat)
    
  }else{stop('check b')}
}
t<-id_fun(b)


clean_exp<-function(data){
  for (i in 1:nrow(data)){
    if ('exposure' %in% names(data)){
      tt<-strsplit(data$exposure[i],split = '||',fixed = TRUE)[[1]][1]
      data$exposure[i]<-trimws(tt, which = c("both", "left", "right"), whitespace = "[ \t\r\n]")
    } else if('outcome' %in% names(data)){
      tt<-strsplit(data$outcome[i],split = '||',fixed = TRUE)[[1]][1]
      data$outcome[i]<-trimws(tt, which = c("both", "left", "right"), whitespace = "[ \t\r\n]")
    } else(print('check columns'))
  }
  return(data)
} 


action<-3

clump_r2 = 0.01
maf_threshold = 0.01

mr_fun<-function(t,id_exposure,exposure_dat){
  for (i in 1:nrow(t$out_list)){
    id_outcome<-t$out_list$id[i]
    print(paste0('MVMR The ',i,' out of ',nrow(t$out_list),' trial!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'))
    print(paste('Exposure ID: ',id_exposure))
    print(paste('Outcome ID: ',id_outcome))
    
    #  folder<-paste(whole_result_folder,'mvmr',id_outcome[i],sep='/')
    
    if (ld_matrix_cal){
      folder<-paste(whole_result_folder,'mvmr_corr',path_name,id_outcome,sep='/')
    } else{
      folder<-paste(whole_result_folder,'mvmr_no_corr',path_name,id_outcome,sep='/')
    }
    dir.create(folder, showWarnings = TRUE, recursive = TRUE)
    setwd(folder)
    save(exposure_dat,file=paste0('exposure_dat',action,'.RData'))
    # if (file.exists(paste0('merge_exp_action',action,'_multi_presso.RData'))){
    #   print(paste0('This ','action ',action,' The ',i,' trial had been done before!'))
    # } else{
    re_log<-paste0('action',action,'_','mvmr_log.log')
    #unlink(re_log,recursive = TRUE)
    if (!develop){
      con <- file(re_log)
      sink(con, append=TRUE) # 记录output
      sink(con, append=TRUE, type="message") # 记录message
    }
    
    
    print(rep('#',20))
    print(paste0('The trial has exp SNP no.: ',nrow(exposure_dat)/2))
    outcome_dat <- extract_outcome_data(exposure_dat$SNP, id_outcome, proxies = TRUE,
                                        maf_threshold = maf_threshold)
    outcome_dat$outcome<-clean_name_fun(outcome_dat$outcome,style=style)
    print(paste0('mr_keep.outcome SNP no.: ', nrow(outcome_dat)))
    table( outcome_dat$mr_keep.outcome)
    outcome_dat<-outcome_dat[outcome_dat$mr_keep.outcome==TRUE,]
    print(paste0('after mr_keep.outcome delete SNP no.: ', nrow(outcome_dat)))
    mvdat <- mv_harmonise_data(exposure_dat, outcome_dat)
    #print(paste0('Ori SNP No.:',nrow(Mydata)))
    save(mvdat,file=paste0('mvmr_har_action',action,'.RData'))
    mv_presso_data<-cbind(mvdat$exposure_beta,mvdat$outcome_beta,mvdat$exposure_se,mvdat$outcome_se)
    mv_presso_data <-as.data.frame(mv_presso_data)
    change_names<-change_name_mvmr(colnames(mv_presso_data),id_outcome)$renam
    out_abbr<-change_name_mvmr(colnames(mv_presso_data),id_outcome)$out_abbr
    colnames(mv_presso_data)<-change_names
    #colnames(mv_presso_data)<-c('RCDW_effect','MCV_effect','Stroke_effect','RCDW_se','MCV_se','Stroke_se')
    
    nrow(distinct(exposure_dat,SNP, .keep_all = TRUE))==nrow(exposure_dat)
    
    data_for_merge<-mv_presso_data
    
    data_for_merge$SNP<-rownames(data_for_merge)
    nrow(distinct(data_for_merge,SNP, .keep_all = TRUE))==nrow(data_for_merge)
    data_for_corr<-data_for_merge%>%inner_join(., distinct(exposure_dat[,c('SNP','effect_allele.exposure','other_allele.exposure')],SNP,.keep_all=TRUE), 
                                               c("SNP" = "SNP"))%>%rename(allele_0=effect_allele.exposure,allele_1=other_allele.exposure)%>%dplyr::select(allele_0,allele_1,RCDW_b,RCDW_se, MCV_b, MCV_se,everything())
    
    stopifnot(nrow(data_for_corr)==nrow(data_for_merge))
    #rownames(data_for_corr)<-data_for_corr$SNP
    
    rawdat_mvmr<-data_for_corr
    save(rawdat_mvmr,file='orginal_har.RData')
    
    
    bx_matrix = as.matrix(dplyr::select(rawdat_mvmr,ends_with('_b'),-paste0(out_abbr,'_b')))
    sebx_matrix = as.matrix(dplyr::select(rawdat_mvmr,ends_with('_se'),-paste0(out_abbr,'_se')))
    #if (file.exists('mr_package_input_snp_correlation_raw.RData'))
    if (ld_matrix_cal){
      if(file.exists('mr_package_input_snp_correlation_raw.RData')){
        load('mr_package_input_snp_correlation_raw.RData')
      } else{
        mr_package_input_snp_corr<-ieugwasr::ld_matrix(
          rawdat_mvmr$SNP,
          plink_bin = genetics.binaRies::get_plink_binary(),
          bfile =bfile)
        save(mr_package_input_snp_corr,file='mr_package_input_snp_correlation_raw.RData')
      }
      
      mr_package_lasso<-mr_mvlasso(mr_mvinput(bx = bx_matrix, bxse = sebx_matrix,
                                              by =  rawdat_mvmr[,paste0(out_abbr,'_b')],correlation = mr_package_input_snp_corr,
                                              byse =rawdat_mvmr[,paste0(out_abbr,'_se')],snps =rawdat_mvmr$SNP,outcome = id_outcome,exposure = id_exposure,effect_allele =rawdat_mvmr$allele_0,other_allele =rawdat_mvmr$allele_1  ))
      
    } else{
      mr_package_lasso<-mr_mvlasso(mr_mvinput(bx = bx_matrix, bxse = sebx_matrix,
                                              by =  rawdat_mvmr[,paste0(out_abbr,'_b')],
                                              byse =  rawdat_mvmr[,paste0(out_abbr,'_se')],snps =rawdat_mvmr$SNP,outcome = id_outcome,exposure = id_exposure,effect_allele =rawdat_mvmr$allele_0,other_allele =rawdat_mvmr$allele_1  ))
      
    }
    
    print(rep(1,20))
    
    ############# #mr_package_input_snp_corr<-ld_matrix(raddat$SNP, with_alleles = TRUE, pop = "EUR")
    
    ####lasso
    rm(mr_package_input_snp_corr)
    gc()
    save(mr_package_lasso,file=paste0('model_action',action,'_MRpackage_lasso.RData'))
    
    print(paste0('Original SNP no.: ', mr_package_lasso@SNPs))
    print(paste0('Valid SNP no.: ', mr_package_lasso@Valid))
    rawdat_mvmr<-rawdat_mvmr[rawdat_mvmr$SNP %in% mr_package_lasso@ValidSNPs,]
    #snp_num<-as.numeric(unlist(strsplit(mr_package_lasso@ValidSNPs,split='snp_',fixed = TRUE)))
    #snp_num<-snp_num[!is.na(snp_num)]
    #rawdat_mvmr<- rawdat_mvmr[snp_num,]
    
    bx_matrix = as.matrix(as.matrix(dplyr::select(rawdat_mvmr,ends_with('_b'),-paste0(out_abbr,'_b'))))
    sebx_matrix = as.matrix(as.matrix(dplyr::select(rawdat_mvmr,ends_with('_se'),-paste0(out_abbr,'_se'))))
    
    save(rawdat_mvmr,file=paste0('Action_',action,'_after_lasso_har_data.RData'))
    # 
    if (ld_matrix_cal){
      if(file.exists('mr_package_input_snp_correlation_after_lasso.RData')){
        load('mr_package_input_snp_correlation_after_lasso.RData')
      }else{
        mr_package_input_snp_corr_after_lasso<-ieugwasr::ld_matrix(
          rawdat_mvmr$SNP,
          plink_bin = genetics.binaRies::get_plink_binary(),
          bfile =bfile
        )
        save(mr_package_input_snp_corr_after_lasso,file='mr_package_input_snp_correlation_after_lasso.RData')
      }
      
      # rmvmr
      mvmr_dat_input<-mr_mvinput(bx = bx_matrix, bxse = sebx_matrix,
                                 by = rawdat_mvmr[,paste0(out_abbr,'_b')],correlation = mr_package_input_snp_corr_after_lasso,
                                 byse =  rawdat_mvmr[,paste0(out_abbr,'_se')],snps =rawdat_mvmr$SNP,outcome = id_outcome,exposure = id_exposure,effect_allele =rawdat_mvmr$allele_0,other_allele =rawdat_mvmr$allele_1  )
      
    } else{
      mvmr_dat_input<-mr_mvinput(bx = bx_matrix, bxse = sebx_matrix,
                                 by = rawdat_mvmr[,paste0(out_abbr,'_b')],
                                 byse = rawdat_mvmr[,paste0(out_abbr,'_se')],snps =rawdat_mvmr$SNP,outcome = id_outcome,exposure = id_exposure,effect_allele =rawdat_mvmr$allele_0,other_allele =rawdat_mvmr$allele_1  )
      
    }
    # rmvmr
    
    # rmvmr_dat <- mrmvinput_to_rmvmr_format(mvmr_dat_input)
    # head(rmvmr_dat)
    # class(rmvmr_dat)
    # 
    # 
    # print(rep('#',20))
    # print(paste0('The ',i,' after har SNP No.:',nrow(rawdat_mvmr)))
    # F.rmvmr<-format_rmvmr(BXGs=dplyr::select(rawdat_mvmr,ends_with('_b'),-paste0(out_abbr,'_b')),
    #                       BYG=rawdat_mvmr[,paste0(out_abbr,'_b')],
    #                       seBXGs=dplyr::select(rawdat_mvmr,ends_with('_se'),-paste0(out_abbr,'_se')),
    #                       seBYG=rawdat_mvmr[,paste0(out_abbr,'_se')],
    #                       RSID=rawdat_mvmr$SNP)
    # if (gencov_cal){
    #   covv<-read.table(corfile,header=F)
    #   covv<-as.matrix(covv)
    #   
    #   phecov<-phenocov_mvmr(covv, rmvmr_dat[,c('sebetaX1' ,  'sebetaX2',  'sebetaX3')])
    #   print(rep('#',10))
    #   print('F')
    #   F_stre_rmvmr<-strength_rmvmr(r_input=rmvmr_dat,gencov=phecov)
    # } else{
    #   print(rep('#',10))
    #   print('F')
    #   F_stre_rmvmr<-strength_rmvmr(F.rmvmr)
    # }
    # 
    # print(F_stre_rmvmr$f)
    # 
    # if (server=='hpc'|server=='featurize'){
    #   pdf(file = "rmvmr_strength_plot_rcdw.pdf", paper = "a4")
    #   print(F_stre_rmvmr$plot[[1]])
    #   dev.off()
    #   pdf(file = "rmvmr_strength_plot_mcv.pdf", paper = "a4")
    #   print(F_stre_rmvmr$plot[[2]])
    #   dev.off()
    # } else {
    #   print(F_stre_rmvmr$plot[[1]])
    #   dev.copy2pdf(file = "rmvmr_strength_plot_rcdw.pdf", paper = "a4r",out.type = "pdf")
    #   dev.off()
    #   print(F_stre_rmvmr$plot[[2]])
    #   dev.copy2pdf(file = "rmvmr_strength_plot_mcv.pdf", paper = "a4r",out.type = "pdf")
    #   dev.off()
    # }
    # 
    # ivw_rmvmr_result<-ivw_rmvmr(r_input = rmvmr_dat, summary = TRUE)
    # save(ivw_rmvmr_result,file=paste0('model_action',action,'_ivw_rmvmr_result.RData'))
    # plot_rmvmr_result<-plot_rmvmr(mvmr_dat_input,ivw_rmvmr_result)
    # if (server=='hpc'|server=='featurize'){
    #   pdf(file = "rmvmr_ivw_plot_wo_correction.pdf", paper = "a4")
    #   print(plot_rmvmr_result$p1)
    #   dev.off()
    #   pdf(file = "rmvmr_ivw_plot_correction.pdf", paper = "a4")
    #   print(plot_rmvmr_result$p2)
    #   dev.off()
    # } else {
    #   print(plot_rmvmr_result$p1)
    #   dev.copy2pdf(file = "rmvmr_ivw_plot_wo_correction.pdf", paper = "a4r",out.type = "pdf")
    #   dev.off()
    #   print(plot_rmvmr_result$p2)
    #   dev.copy2pdf(file = "rmvmr_ivw_plot_correction.pdf", paper = "a4r",out.type = "pdf")
    #   dev.off()
    # }
    # pleiotropy_rmvmr_result<-  pleiotropy_rmvmr(mvmr_dat_input,ivw_rmvmr_result)
    #####other model
    
    
    ####ivw fe
    mvmr_ivw_fe = mr_mvivw(mvmr_dat_input,model='fixed')
    save(mvmr_ivw_fe,file=paste0('model_action',action,'_MRpackage_ivw_fe.RData'))
    ####ivw re
    mvmr_ivw_re = mr_mvivw(mvmr_dat_input,model='random')
    save( mvmr_ivw_re,file=paste0('model_action',action,'_MRpackage_ivw_re.RData'))
    ### median
    set.seed(20200706)
    mvmr_MR_median<-mr_mvmedian(mvmr_dat_input,
                                seed = 3
    )
    save(mvmr_MR_median,file=paste0('model_action',action,'_MRpackage_median.RData'))
    # egger
    mvmr_egger_AD = mr_mvegger(mvmr_dat_input)
    save(mvmr_egger_AD,file=paste0('model_action',action,'_egger.RData'))
    
    # robust
    mvmr_robust_AD = mvmr_robust(bx_matrix, rawdat_mvmr[,paste0(out_abbr,'_b')],
                                 rawdat_mvmr[,paste0(out_abbr,'_se')], k.max = 1000, maxit.scale = 1000)
    save(mvmr_robust_AD,file=paste0('model_action',action,'_robust.RData'))
    dat = data.frame("by" = rawdat_mvmr[,paste0(out_abbr,'_b')], "seby" =rawdat_mvmr[,paste0(out_abbr,'_se')],
                     "bx" = bx_matrix, "sebx" = sebx_matrix)
    
    #  median
    set.seed(20200706)
    mvmr_median_AD = mvmr_median(bx_matrix, sebx_matrix, rawdat_mvmr[,paste0(out_abbr,'_b')],
                                 rawdat_mvmr[,paste0(out_abbr,'_se')], boot = TRUE, boot_it = 1000)
    save(mvmr_median_AD,file=paste0('model_action',action,'_median.RData'))
    #自制模型 lasso
    set.seed(20200707)
    mvmr_lasso_AD = mvmr_lasso(bx_matrix,rawdat_mvmr[,paste0(out_abbr,'_b')],
                               rawdat_mvmr[,paste0(out_abbr,'_se')])
    
    save(mvmr_lasso_AD,file=paste0('model_action',action,'_lasso.RData'))
    
    
    
    
    
    
    # go on to build other models
    
    print(rep('#',20))
    print(paste0('The ',i,' after har SNP No.:',nrow(rawdat_mvmr)))
    F.data<-format_mvmr(BXGs=dplyr::select(rawdat_mvmr,ends_with('_b'),-paste0(out_abbr,'_b')),
                        BYG=rawdat_mvmr[,paste0(out_abbr,'_b')],
                        seBXGs=dplyr::select(rawdat_mvmr,ends_with('_se'),-paste0(out_abbr,'_se')),
                        seBYG=rawdat_mvmr[,paste0(out_abbr,'_se')],
                        RSID=rawdat_mvmr$SNP)
    
    
    
    save(F.data,file=paste0('action',action,'_har','.RData'))
    head(F.data)
    if (gencov_cal){
      covv<-read.table(corfile,header=F)
      covv<-as.matrix(covv)
      
      phecov<-phenocov_mvmr(covv, F.data[,c('sebetaX1' ,  'sebetaX2',  'sebetaX3')])
      print(rep('#',10))
      print('F')
      F_stre<-strength_mvmr(r_input=F.data,gencov=phecov)
      print(rep('#',10))
      #  print('F by minimisation of Q-statistics')
      #  F_qhe<-strhet_mvmr(r_input=F.data,gencov=phecov)
      print(rep('#',10))
      print('Pleiotropy')
      plei<-pleiotropy_mvmr(r_input=F.data,gencov=phecov)
      print(rep('#',10))
      print('ivw mvmr')
      mvmr_ivw<-ivw_mvmr(r_input=F.data,gencov = covv)
      save(mvmr_ivw,file=paste0('merge_exp_action',action,'_mvmr_ivw.RData'))
      print(rep('#',10))
    }else{
      print(rep('#',10))
      print('F')
      F_stre<-strength_mvmr(r_input=F.data,gencov=0)
      print(rep('#',10))
      #  print('F by minimisation of Q-statistics')
      #  F_qhe<-strhet_mvmr(r_input=F.data,gencov=phecov)
      print(rep('#',10))
      print('Pleiotropy')
      plei<-pleiotropy_mvmr(r_input=F.data,gencov=0)
      print(rep('#',10))
      print('ivw mvmr')
      mvmr_ivw<-ivw_mvmr(r_input=F.data,gencov = 0)
      save(mvmr_ivw,file=paste0('merge_exp_action',action,'_mvmr_ivw.RData'))
      print(rep('#',10))
    }
    vars<-gsub('_b','',names(dplyr::select(rawdat_mvmr,ends_with('_b'))))
    F_stats=list('F_stre'=F_stre,'plei'=plei,'vars'=vars)
    #  F_stats=list('F_stre'=F_stre,'F_qhe'=F_qhe,'plei'=plei,'F_stre_rmvmr'=F_stre_rmvmr,'pleiotropy_rmvmr_result'=pleiotropy_rmvmr_result)
    save(F_stats,file=paste0('action',action,'_F_stats.RData'))
    
    
    
    
    print(rep('#',10))
    print('TSMR mvmr')
    #pdf(file = "mv_mutiple.pdf")
    res_mv_mutiple <- mv_multiple(mvdat,plots = TRUE)
    res_mv_mutiple_intecept<-mv_multiple(
      mvdat,
      intercept = TRUE,
      plots = TRUE
    )
    #print(res_mv_mutiple)
    # dev.off()
    
    save(res_mv_mutiple ,file=paste0('merge_exp_action',action,'_res_mv_mutiple.RData'))
    save(res_mv_mutiple_intecept ,file=paste0('merge_exp_action',action,'_res_mv_mutiple_intecept.RData'))
    if (server=='hpc'|server=='featurize'){
      pdf(file = "mv_mutiple1.pdf", paper = "a4")
      print(res_mv_mutiple$plots[[1]])
      dev.off()
      pdf(file = "mv_mutiple2.pdf", paper = "a4")
      print(res_mv_mutiple$plots[[2]])
      dev.off()
      pdf(file = "mv_mutiple_intecept1.pdf", paper = "a4")
      print(res_mv_mutiple_intecept$plots[[1]])
      dev.off()
      pdf(file = "mv_mutiple_intecept2.pdf", paper = "a4")
      print(res_mv_mutiple_intecept$plots[[2]])
      dev.off()
    } else {
      print(res_mv_mutiple$plots[[1]])
      dev.copy2pdf(file = "mv_mutiple1.pdf", paper = "a4r",out.type = "pdf")
      dev.off()
      print(res_mv_mutiple$plots[[2]])
      dev.copy2pdf(file = "mv_multiple2.pdf", paper = "a4r",out.type = "pdf")
      dev.off()
      print(res_mv_mutiple_intecept$plots[[1]])
      dev.copy2pdf(file = "mv_mutiple_intecept1.pdf", paper = "a4r",out.type = "pdf")
      dev.off()
      print(res_mv_mutiple_intecept$plots[[2]])
      dev.copy2pdf(file = "mv_multiple_intecept2.pdf", paper = "a4r",out.type = "pdf")
      dev.off()
    }
    
    
    print(rep('#',10))
    print('presso')
    data_presso<-rawdat_mvmr
    rownames(data_presso)<-NULL
    set.seed(20200705)
    file_presso<-paste0('merge_exp_action',action,'_multi_presso_nb_5000.RData')
    if(file.exists(file_presso)){
      load(file_presso)
    }else{
      presso<-mr_presso(BetaOutcome =paste0(out_abbr,'_b'), BetaExposure = names(dplyr::select(rawdat_mvmr,ends_with('_b'),-paste0(out_abbr,'_b'))), SdOutcome =paste0(out_abbr,'_se'),
                        SdExposure = names(dplyr::select(rawdat_mvmr,ends_with('_se'),-paste0(out_abbr,'_se'))), OUTLIERtest = TRUE,DISTORTIONtest = TRUE,
                        data =data_presso, NbDistribution = 5000, SignifThreshold = 0.05)
      save(presso,file=file_presso)
    }
    
    print(presso)
    if (!develop){
      sink()
      sink(type="message")
      cat(readLines(re_log), sep="\n")
    }
  }
}

if (b=='b2'){
  for (j in 1:nrow(t$med_list)){
    id_exposure<-c(t$exp_id,t$med_list$id[j])
    path_name<-paste(c(t$exp_list$abbr,t$med_list$abbr[j]),collapse='_')
    print(id_exposure)
    print('\n')
    exposure_dat <- mv_extract_exposures(id_exposure,harmonise_strictness =3,clump_r2 = clump_r2)
    nrow(distinct(exposure_dat,SNP, .keep_all = TRUE))==nrow(exposure_dat)
    exposure_dat$exposure<-clean_name_fun(exposure_dat$exposure,style=style)
    
    mr_fun(t,id_exposure,exposure_dat)
  }
}else if((b=='b0')|(b=='b1')){
  id_exposure<-t$exp_id
  print(id_exposure)
  print('\n')
  exposure_dat <- mv_extract_exposures(id_exposure,harmonise_strictness =3,clump_r2 = clump_r2)
  nrow(distinct(exposure_dat,SNP, .keep_all = TRUE))==nrow(exposure_dat)
  exposure_dat$exposure<-clean_name_fun(exposure_dat$exposure,style=style)
  path_name<-paste(t$exp_list$abbr,collapse='_')
  mr_fun(t,id_exposure,exposure_dat)
}





###################################################################################################



####################################################################################################################


print(paste0('time spent:',' ',proc.time()-ptm))
print('Done!!!!')
