
library("readxl")
library(TwoSampleMR)
library(phenoscanner)
library(LDlinkR)
library(dplyr)
library(cause)
library(RadialMR)
library(tidyverse)
library(data.table)
library(MendelianRandomization)
#dev.off()
rm(list=ls())
gc()

server='home'
whole_result_folder<-'.../2_20/result/umr'
data_table<-'.../2_20/data_table.csv'
data_list<-list.files('.../2_20/result/umr',pattern='_har.RData$',full.names=TRUE,recursive=TRUE)
source('.../2_20/code/umr_tool.R')

data_list<-data.frame(data_list)
data_list$folder<-dirname(data_list$data_list)

path_fun0<-function(data){
  re<-data[length(data)]
  return(re)
}
path_fun1<-function(data){
  re<-data[length(data)-1]
  return(re)
}
path_fun2<-function(data){
  re<-data[length(data)-2]
  return(re)
}
data_list$out_id<-unlist(lapply(strsplit(data_list$folder,split='/'),path_fun0))
data_list$exp_id<-unlist(lapply(strsplit(data_list$folder,split='/'),path_fun1))
data_list<-as_tibble(data_list)
data_table<-read_csv(data_table)
final<-inner_join(data_list, data_table[,c('id','trait')], 
           c("out_id" = "id"))%>%rename(out_trait = trait)%>%inner_join(., data_table[,c('id','trait')], 
                                                                     c("exp_id" = "id"))%>%rename(exp_trait = trait)
final<-final[!endsWith(final$data_list, 'final_har.RData'),]
ptm <- proc.time()
whole_re<-data.frame()
whole_re_forest<-data.frame()


or_fun<-function (data1,all_res_forest,folder,raddat,type='presso') 
{
  if (type=='presso'){
    mr_res<-all_res_forest[1:2,]
    mr_res$Method<-c('MR PRESSO','MR PRESSO(outlier corrected)')
    mr_res$b<-data1$`Main MR results`$`Causal Estimate`
    mr_res$se<-data1$`Main MR results`$Sd
    mr_res$pval<-data1$`Main MR results`$`P-value`
    mr_res$lo_ci <- mr_res$b - 1.96 * mr_res$se
    mr_res$up_ci <- mr_res$b + 1.96 * mr_res$se
    mr_res$or <- exp(mr_res$b)
    mr_res$or_lci95 <- exp(mr_res$lo_ci)
    mr_res$or_uci95 <- exp(mr_res$up_ci)
    mr_res$lo_ci<-NULL
    mr_res$up_ci<-NULL
    mr_res[,c("Q", "Q_df", "Q_pval", "intercept", "intercept_se", "intercept_pval",'SNP')]<-NA
    mr_res$presso_pval<-c(data1$`MR-PRESSO results`$`Global Test`$Pvalue,NA)
  } else  if (type=='rad_ivw'){
    
    ivw_log_folder<-paste(folder,'rad_ivw.txt',sep='/')
    unlink(ivw_log_folder,recursive = TRUE)
    con_ivw <- file(ivw_log_folder) # 创建一???.log文件
    sink(con_ivw, append=TRUE) # 记录output
    sink(con_ivw, append=TRUE, type="message") # 记录message
    ivw_radial(raddat,0.05/nrow(raddat),3)
    sink()
    sink(type="message")
    cat(readLines(ivw_log_folder), sep="\n")
    ivw_txt<-read.delim(ivw_log_folder,sep = "\t")
    unlink(ivw_log_folder,recursive = TRUE)
    #    ivw_Q_p<-as.numeric(strsplit(a$Radial.IVW[8],split='p-value: ', fixed = TRUE)[[1]][2])
    ivw_Q<-as.numeric(str_extract_all(ivw_txt$Radial.IVW[8],"\\d+\\.\\d+|\\d+")[[1]])[1]
    ivw_Q_df<-as.numeric(str_extract_all(ivw_txt$Radial.IVW[8],"\\d+\\.\\d+|\\d+")[[1]])[2]
    ivw_Q_p<-as.numeric(str_extract_all(ivw_txt$Radial.IVW[8],"\\d+\\.\\d+|\\d+")[[1]])[3]
    
    mr_res<-all_res_forest[1:2,]
    mr_res$Method<-c('RadialMR IVW(FE)','RadialMR IVW(RE)')
    mr_res$b<-data1$coef$Estimate[3:4]
    mr_res$se<-data1$coef$Std.Error[3:4]
    mr_res$pval<-data1$coef$`Pr(>|t|)`[3:4]
    mr_res$lo_ci <- mr_res$b - 1.96 * mr_res$se
    mr_res$up_ci <- mr_res$b + 1.96 * mr_res$se
    mr_res$or <- exp(mr_res$b)
    mr_res$or_lci95 <- exp(mr_res$lo_ci)
    mr_res$or_uci95 <- exp(mr_res$up_ci)
    mr_res$lo_ci<-NULL
    mr_res$up_ci<-NULL
    mr_res$nsnp<-c(length(raddat$SNP),length(raddat$SNP))
    mr_res$Q<-c(ivw_Q,ivw_Q)
    mr_res$Q_df<-c(ivw_Q_df,ivw_Q_df)
    mr_res$Q_pval<-c(ivw_Q_p,ivw_Q_p)
    mr_res[,c("intercept", "intercept_se", "intercept_pval",'SNP','presso_pval')]<-NA
  }
  else  if (type=='rad_egger'){
    egger_log_folder<-paste(folder,'rad_egger.txt',sep='/')
    unlink(egger_log_folder,recursive = TRUE)
    con_egger <- file(egger_log_folder) # 创建一???.log文件
    sink(con_egger, append=TRUE) # 记录output
    sink(con_egger, append=TRUE, type="message") # 记录message
    egger_radial(raddat,0.05/nrow(raddat),3)
    sink()
    sink(type="message")
    cat(readLines(egger_log_folder), sep="\n")
    egger_txt<-read.delim(egger_log_folder,sep = "\t")
    unlink(egger_log_folder,recursive = TRUE)
    egger_Q<-as.numeric(str_extract_all(egger_txt$Radial.MR.Egger[6],"\\d+\\.\\d+|\\d+")[[1]])[1]
    egger_Q_df<-as.numeric(str_extract_all(egger_txt$Radial.MR.Egger[6],"\\d+\\.\\d+|\\d+")[[1]])[2]
    egger_Q_p<-as.numeric(str_extract_all(egger_txt$Radial.MR.Egger[6],"\\d+\\.\\d+|\\d+")[[1]])[3]
    egger<-data.frame(data1$coef)
    names(egger)<-c('b','se','t','p')
    
    mr_res<-all_res_forest[1,]
    mr_res$Method<-c('RadialMR Egger')
    mr_res$b<-egger$b[2]
    mr_res$se<-egger$se[2]
    mr_res$pval<-egger$p[2]
    mr_res$lo_ci <- mr_res$b - 1.96 * mr_res$se
    mr_res$up_ci <- mr_res$b + 1.96 * mr_res$se
    mr_res$or <- exp(mr_res$b)
    mr_res$or_lci95 <- exp(mr_res$lo_ci)
    mr_res$or_uci95 <- exp(mr_res$up_ci)
    mr_res$lo_ci<-NULL
    mr_res$up_ci<-NULL
    mr_res$nsnp<-c(length(raddat$SNP))
    mr_res$Q<-egger_Q
    mr_res$Q_df<-egger_Q_df
    mr_res$Q_pval<-egger_Q_p
    mr_res$intercept_pval<-egger$p[1]
    mr_res$intercept<-egger$b[1]
    mr_res$intercept_se<-egger$se[1]
    
    mr_res[,c('SNP')]<-NA
    mr_res$presso_pval<-NA
  } else if (type=='mr_package_ivw'){
    mr_res<-all_res_forest[1:2,]
    mr_res$Method<-c('MR IVW(FE)','MR IVW(RE)')
    mr_res$nsnp<-c(mr_package_ivw_fe@SNPs,mr_package_ivw_re@SNPs)
    mr_res$b<-c(mr_package_ivw_fe@Estimate,mr_package_ivw_re@Estimate)

    mr_res$se<-c(mr_package_ivw_fe@StdError,mr_package_ivw_re@StdError)
    mr_res$pval<-c(mr_package_ivw_fe@Pvalue,mr_package_ivw_re@Pvalue)
    mr_res$lo_ci <-c(mr_package_ivw_fe@CILower,mr_package_ivw_re@CILower)
    mr_res$up_ci <-c(mr_package_ivw_fe@CIUpper,mr_package_ivw_re@CIUpper)
    mr_res$or <- exp(mr_res$b)
    mr_res$or_lci95 <- exp(mr_res$lo_ci)
    mr_res$or_uci95 <- exp(mr_res$up_ci)
    mr_res$lo_ci<-NULL
    mr_res$up_ci<-NULL
    mr_res$nsnp<-c(mr_package_ivw_fe@SNPs,mr_package_ivw_re@SNPs)
    mr_res$Q<-c(mr_package_ivw_fe@Heter.Stat[1],mr_package_ivw_re@Heter.Stat[1])
    mr_res$Q_pval<-c(mr_package_ivw_fe@Heter.Stat[2],mr_package_ivw_re@Heter.Stat[2])
    mr_res[,c( "Q_df", "intercept", "intercept_pval","intercept_se", 'SNP','lower_intercept', 'upper_intercept')]<-NA
    
  } else if (type=='mr_package_egger'){
    mr_res<-all_res_forest[1,]
    mr_res$Method<-c('MR Egger')
    mr_res$nsnp<-c(data1@SNPs)
    mr_res$b<-c(data1@Estimate)
    
    mr_res$se<-c(data1@StdError.Est)
    mr_res$pval<-c(data1@Pvalue.Est)
    mr_res$lo_ci <-c(data1@CILower.Est)
    mr_res$up_ci <-c(data1@CIUpper.Est)
    mr_res$or <- exp(mr_res$b)
    mr_res$or_lci95 <- exp(mr_res$lo_ci)
    mr_res$or_uci95 <- exp(mr_res$up_ci)
    mr_res$lo_ci<-NULL
    mr_res$up_ci<-NULL
    mr_res$nsnp<-c(data1@SNPs)
    mr_res$Q<-c(data1@Heter.Stat[1])
    mr_res$Q_pval<-c(data1@Heter.Stat[2])
    mr_res$intercept<-c(data1@Intercept)
    mr_res$intercept_se<-c(data1@StdError.Int)
    mr_res$intercept_pval<-c(data1@Pvalue.Int)
    mr_res$lower_intercept <-c(data1@CILower.Int)
    mr_res$upper_intercept <-c(data1@CIUpper.Int)
    mr_res[,c( "Q_df",'SNP')]<-NA
  }
  
  return(mr_res)
}

# i<-1
# l<-list()
# l_other<-list()
for (i in 10:12){#1:length(final$data_list)
  
  print(paste('The ', i, ' trial'))
  folder<-final$folder[i]
  setwd(folder)
#  re_log_folder<-'radial_outliner.log'
#  unlink(re_log_folder,recursive = TRUE)
  #  dir.create(f, showWarnings = TRUE, recursive = TRUE)
#  con <- file(re_log_folder) # 
#  sink(con, append=TRUE) # 
#  sink(con, append=TRUE, type="message") # 
  load(final$data_list[i])
  exp_id<-final$exp_id[i]
  out_id<-final$out_id[i]
  exp_trait<-final$exp_trait[i]
  out_trait<-final$out_trait[i]
  print(paste('Outcome: ',out_trait,'; ID ',out_id))
  print(paste0('Exp: ', har$exposure[1]))
  print(paste0('Out: ',har$outcome[1]))
  print(paste0('Ori SNP No.:',nrow(har)))
  dat<-har
  dat<-dat[dat["ambiguous" ]==FALSE,]
  print(paste0('Ori SNP No.after deleting ambiguous:',nrow(dat)))
  
  raddat <- format_radial(BXG=dat$beta.exposure, BYG=dat$beta.outcome, 
                          seBXG=dat$se.exposure, seBYG=dat$se.outcome, RSID=dat$SNP)
  outliers<-'NA'
  outliner_trial<-0
  while (outliers!="No significant outliers") {
    outliner_trial=outliner_trial+1
    ivwrad <- ivw_radial(raddat, alpha=0.05, weights=3)
    outliers=ivwrad$outliers[1]
    if (outliers!="No significant outliers"){
      myvars=outliers$SNP
      print(paste0('Radial IVW ',outliner_trial,' trail has ',length(myvars), ' outliners without p adjusting: ',myvars))
      raddat <- raddat[ ! raddat$SNP %in% myvars, ]
      print(rep('1',15))
    } else {
      print(paste0('Radial IVW ',outliner_trial,' trail has NOT outliners without p adjusting'))
      raddat<-raddat
      print(rep('1——end',30))
    }
  }
  
  outliers2<-'NA'
  outliner_trial2<-0
  while (outliers2!="No significant outliers") {
    
    outliner_trial2=outliner_trial2+1
    egger2 <- egger_radial(raddat,0.05,3)
    outliers2=egger2 $outliers[1]
    if (outliers2!="No significant outliers"){
      myvars=outliers2$SNP
      raddat <- raddat[ ! raddat$SNP %in% myvars, ]
      print(rep('2',15))
      print(paste0('Radial Egger ',outliner_trial,' trail has',length(myvars), ' outliners after p adjusting: ',myvars))
    } else {
      print(paste0('Radial Egger ',outliner_trial,' trail has NOT outliners after p adjusting'))
      raddat<-raddat
      print(rep('2--',30))
    }
  }
  
  
  
  
  
  dat2 <- dat[  dat$SNP %in% raddat$SNP, ]
  print(paste0('Before Radial has ',nrow(dat),' snps; Then has ', nrow(dat2), ' snps.'))
  dat2<-steiger_filtering(dat2)#%>%subset(steiger_dir==TRUE)
  print(paste0('After Steiger has ', nrow(dat2), ' snps.'))
  
  
  
  
  singlesnp_results=mr_singlesnp(dat2, parameters = default_parameters(),
                                 single_method = "mr_wald_ratio", 
                                 all_method = c("mr_ivw",
                                                "mr_egger_regression", 
                                                "mr_weighted_mode",
                                                "mr_weighted_median"))
  
  single_snp_all<-singlesnp_results[singlesnp_results$SNP %in% grep("^All", singlesnp_results$SNP, value = T),]
  single_snp_single<-singlesnp_results[singlesnp_results$SNP %in% grep("^rs", singlesnp_results$SNP, value = T),]
  if (all(single_snp_all$p>=0.05)&any(single_snp_single$p<0.05)){
    print('Outliners in single snp mr result:')
    print(singlesnp_results[single_snp_single$p<0.05,][,c('SNP','p')])
    out_snp<-singlesnp_results[single_snp_single$p<0.05,]$SNP
    dat2<-dat2[which(!(dat2$SNP %in% out_snp)),]
    
  } else if (all(single_snp_all$p<0.05)&any(single_snp_single$p>=0.05)){
    print('Outliners in single snp mr result:')
    print(singlesnp_results[single_snp_single$p>=0.05,][,c('SNP','p')])
    out_snp<-singlesnp_results[single_snp_single$p>=0.05,]$SNP
    dat2<-dat2[which(!(dat2$SNP %in% out_snp)),]
  }
  print(paste0('After single sno test has ', nrow(dat2), ' snps.'))
  rownames(dat2)<-NULL
  
  raddat <- raddat[   raddat$SNP %in% dat2$SNP, ]
  
  save(raddat,file=paste(folder,'final_har.RData',sep='/'))
  
  print('Final Radial IVW')
  ivwrad3 <- ivw_radial(raddat,0.05/nrow(raddat),3)
  print('Final Radial egger')
  egger.model3<-egger_radial(raddat,0.05/nrow(raddat),3)
  #print(egger.model3)
  # mr package
  print(rep(0,20))
  if (file.exists('mr_package_input_snp_correlation.RData')){
    load('mr_package_input_snp_correlation.RData')
  }else{
    mr_package_input_snp_corr<-ieugwasr::ld_matrix(
      raddat$SNP,
      plink_bin = genetics.binaRies::get_plink_binary(),
      bfile = ".../plink_reference/EUR"
    )
    save(mr_package_input_snp_corr,file='mr_package_input_snp_correlation.RData')
    
  }
   print(rep(1,20))
  
  #mr_package_input_snp_corr<-ld_matrix(raddat$SNP, with_alleles = TRUE, pop = "EUR")
  mr_package_input<-mr_input(bx = raddat$beta.exposure, bxse = raddat$se.exposure, by = raddat$beta.outcome, byse = raddat$se.outcome,snps = raddat$SNP,correlation=mr_package_input_snp_corr)
 # # rm(mr_package_input_snp_corr)
 #  gc()
  mr_package_lasso<-mr_lasso(mr_package_input)
  if (mr_package_lasso@SNPs!=mr_package_lasso@Valid){
    print(paste('Outcome: ',out_trait,'; ID ',out_id))
    print(paste0('Exp: ', har$exposure[1]))
    print(paste0('Out: ',har$outcome[1]))
    stop('Lasso SNPS have outliers!')
  }
  # print(rep(2,20))   
  # mr_package_ivw_re<-mr_ivw(mr_package_input,correl = TRUE,model ='random')
  # mr_package_ivw_fe<-mr_ivw(mr_package_input,correl = TRUE,model ='fixed')
  # print(rep(3,20))
  # mr_package_egger<-mr_egger(mr_package_input,correl=TRUE)
  #######################################
  if (server=='hpc'|server=='featurize'){
    pdf(file = "IVW_T_T_F.pdf", paper = "a4")
    IVWplot1<-plot_radial(ivwrad3,T,T,F)
    print(IVWplot1)
    dev.off()
    
    pdf(file = "IVW_T_T_T.pdf", paper = "a4")
    IVWplot2<-plot_radial(ivwrad3,T,T,T)
    print(IVWplot2)
    dev.off()

  } else {
    IVWplot1<-plot_radial(ivwrad3,T,T,F)
    print(IVWplot1)
    dev.copy2pdf(file = "IVW_T_T_F.pdf", paper = "a4",out.type = "pdf")
    dev.off()
    IVWplot2<-plot_radial(ivwrad3,T,T,T)
    print(IVWplot2)
    dev.copy2pdf(file = "IVW_T_T_T.pdf", paper = "a4",out.type = "pdf")
    dev.off()
  }
  
  if (server=='hpc'|server=='featurize'){
    egger.model<-egger_radial(raddat,0.05/nrow(raddat),3)
    pdf(file = "eggerT_T_F.pdf", paper = "a4")
    Eggerplot1<-plot_radial(egger.model,T,T,F)
    print(Eggerplot1)
    dev.off()
    
    pdf(file = "egger_T_T_T.pdf", paper = "a4")
    Eggerplot2<-plot_radial(egger.model,T,T,T)
    print(Eggerplot2)
    dev.off()
  } else {
    egger.model<-egger_radial(raddat,0.05/nrow(raddat),3)
    Eggerplot1<-plot_radial(egger.model,T,T,F)
    print(Eggerplot1)
    dev.copy2pdf(file = "eggerT_T_F.pdf", paper = "a4",out.type = "pdf")
    dev.off()
    Eggerplot2<-plot_radial(egger.model,T,T,T)
    print(Eggerplot2)
    dev.copy2pdf(file = "egger_T_T_T.pdf", paper = "a4",out.type = "pdf")
    dev.off()
    }

  
  #run MR without outliers
  mr_results <- mr(dat2)
  res<-mr_results
  #l[[i]]<-res
  print(rep('#',10))
  print('MR result')
  print(res)
  
  het<-mr_heterogeneity(dat2, method_list=c("mr_egger_regression", "mr_ivw"))  #
  print(rep('#',10))
  print('Heterogeneity test')  
  print(het)
  # 
  plt <- mr_pleiotropy_test(dat2) ####egger
  print(rep('#',10))
  print('Pleiotropy test')  
  print(plt)
  
  
  
  
  library(MRPRESSO)
  presso<-mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data =dat2, NbDistribution = 1000, SignifThreshold = 0.05)
  print(rep('#',10))
  print('presso')
  print(presso)

  

  if (server=='hpc'|server=='featurize'){
    single <- mr_leaveoneout(dat2)
    single
    pdf(file = "leave_one_out.pdf", paper = "a4")
    print(mr_leaveoneout_plot(single))   #
    dev.off()

    pdf(file = "scatter.pdf", paper = "a4")
    print(mr_scatter_plot(res,dat2))#
    dev.off()

    pdf(file = "forest.pdf", paper = "a4")
    res_single <- mr_singlesnp(dat2)
    print(mr_forest_plot(res_single))#
    dev.off()

    pdf(file = "funnel.pdf", paper = "a4")
    print(mr_funnel_plot(res_single))#
    dev.off()

    res_odd<-generate_odds_ratios(res)#
  } else {
    single <- mr_leaveoneout(dat2)
    single
    print(mr_leaveoneout_plot(single))   #
    dev.copy2pdf(file = "leave_one_out.pdf", paper = "a4",out.type = "pdf")
    dev.off()
    print(mr_scatter_plot(res,dat2))#
    dev.copy2pdf(file = "scatter.pdf", paper = "a4",out.type = "pdf")
    dev.off()
    res_single <- mr_singlesnp(dat2)
    print(mr_forest_plot(res_single))#
    dev.copy2pdf(file = "forest.pdf", paper = "a4",out.type = "pdf")
    dev.off()
    print(mr_funnel_plot(res_single))#
    dev.copy2pdf(file = "funnel.pdf", paper = "a4",out.type = "pdf")
    dev.off()
    res_odd<-generate_odds_ratios(res)#
  }

  
  singlesnp_results=mr_singlesnp(dat2, parameters = default_parameters(),
                                 single_method = "mr_wald_ratio", 
                                 all_method = c("mr_ivw",
                                                "mr_egger_regression", 
                                                "mr_weighted_mode",
                                                "mr_weighted_median"))
  print(rep('#',10))
  singlesnp_results.2=grep("^rs", singlesnp_results$SNP, value = T)
  print(paste0('Single SNP method No.: ',length(singlesnp_results.2)))
  dat3=dat2[which(dat2$SNP %in% singlesnp_results.2),]
  dim(dat3)
  

  Q_ivw<-het$Q[1]
  Q_ivw_df<-het$Q_df[1]
  print(paste('Corhan Q I2 ',round((Q_ivw-Q_ivw_df)/Q_ivw,4)))
  
  
  print(rep('#',10))
  print('F test') 
  #r2 and Fstat
  #removed 2 NA's that didn't have EAF
  n=dat2$samplesize.exposure[1]
  k=nrow(dat2)
  p=dat2$eaf.exposure
  r2.1<-(2*(dat2$beta.exposure^2)*p*(1-p))/1
  r2=sum(r2.1, na.rm = TRUE)
  print(paste0('R^2: ',r2))
  Fstat <- r2*(n-1-k)/((1-r2)*k)
  print(paste0('F: ',Fstat))
  
  maf_I2<-round((Q_ivw-Q_ivw_df)/Q_ivw,4)
  
  print(rep('#',20))
  print('combined results') 
  #Combine results
  sin=singlesnp_results
  res=mr_results
  het<-mr_heterogeneity(dat2)
  plt<-mr_pleiotropy_test(dat2)
  
  all_res<-nn(res,het,plt,sin,ao_slc=T,
              Exp=T,split.exposure=T,split.outcome=T)
  all_res$lower_intercept=all_res$intercept-1.96*all_res$intercept_se
  all_res$upper_intercept=all_res$intercept+1.96*all_res$intercept_se
  all_res$presso_pval<-NA
  
  all_res_forest<-subset(all_res,Method!='Wald ratio')
  rad_ivw<-or_fun(ivwrad3,all_res_forest,folder,raddat,type='rad_ivw')
  rad_egger<-or_fun(egger.model3,all_res_forest,folder,raddat,type='rad_egger')
  presso_re<-or_fun(presso,all_res_forest,folder,raddat,type='presso')
  print(rep('X',50))
  all_res<-rbind(all_res,rad_ivw,rad_egger,presso_re)#,mr_package_ivw_result,mr_package_egger_result)
  
  # mr_package_ivw_result<-or_fun(NA,all_res_forest,folder,raddat,type='mr_package_ivw')
  # mr_package_egger_result<-or_fun(mr_package_egger,all_res_forest,folder,raddat,type='mr_package_egger')


  all_res_forest1<-rbind(all_res_forest,rad_ivw,rad_egger,presso_re)#,mr_package_ivw_result,mr_package_egger_result)
  
  all_res_forest1$r2_maf_sum<-r2_maf_sum
  all_res_forest1$r2_maf_sum[(length(all_res_forest1$r2_maf_sum)-1):length(all_res_forest1$r2_maf_sum)]<-NA
  all_res_forest1$F_maf<-Fstat_maf
  all_res_forest1$F_maf[(length(all_res_forest1$F_maf)-1):length(all_res_forest1$F_maf)]<-NA
  all_res_forest1$maf_I2<-maf_I2
  all_res_forest1$maf_I2[(length(all_res_forest1$maf_I2)-1):length(all_res_forest1$maf_I2)]<-NA
  all_res_forest1$r2_eaf<-r2
  all_res_forest1$r2_eaf[(length(all_res_forest1$r2_eaf)-1):length(all_res_forest1$r2_eaf)]<-NA
  all_res_forest1$F_eaf<-Fstat
  all_res_forest1$F_eaf[(length(all_res_forest1$F_eaf)-1):length(all_res_forest1$F_eaf)]<-NA

  # 
  # ivw_fe_logic<-tryCatch({
  #   mr_package_ivw_fe<-mr_ivw(mr_package_input,correl = TRUE,model ='fixed')
  # }, warning = function(w){
  #   print(FALSE)
  # }, error = function(e){
  #   print(FALSE)
  # },finally = {
  #   print(TRUE)
  # })
  # ivw_re_logic<-tryCatch({
  #   mr_ivw(mr_package_input,correl = TRUE,model ='random')
  # }, warning = function(w){
  #   print(FALSE)
  # }, error = function(e){
  #   print(FALSE)
  # },finally = {
  #   print(TRUE)
  # })
  # 
  # if (class(ivw_fe_logic)!="logical"){
  #   print(rep('Y',50))
  #   model_result<-mr_package_ivw_fe
  #   fe_log_folder<-paste(folder,'mr_package_ivw.Rout',sep='/')
  #   unlink(fe_log_folder,recursive = TRUE)
  #   con_fe <- file(fe_log_folder, open = "wt") # 
  #   sink(con_fe) # 
  #   sink(con_fe, append=TRUE, type="message") # 
  #   print(try(mr_ivw(mr_package_input,correl = TRUE,model ='random'),silent = TRUE))
  #   sink()
  #   sink(type="message")
  #   file.show(fe_log_folder)
  #   mr_fe_I2<-read.delim(fe_log_folder,sep = '\t')
  #   unlink(fe_log_folder,recursive = TRUE)
  #   mr_fe_I2<-mr_fe_I2$Inverse.variance.weighted.method[nrow(mr_fe_I2)]
  #   mr_fe_I2<-gsub("[[:space:]]", "", strsplit(mr_fe_I2,split='=')[[1]][4])
  #   mr_fe_I2<-as.numeric(substr(mr_fe_I2,1,nchar(mr_fe_I2)-2))/100
  #   if (is.na(mr_fe_I2)){
  #     mr_fe_I2<-'NA'
  #   }
  #   all_res_forest1[all_res_forest1$Method=="MR IVW(FE)",'maf_I2']<-mr_fe_I2
  #   all_res_forest1[all_res_forest1$Method=="MR IVW(FE)",'r2_eaf']<-mr_fe_I2
  # }
  # 
  # if (class(ivw_re_logic)!="logical"){
  # print(rep('z',50))
  # model_result<-mr_package_ivw_re
  # fe_log_folder<-paste(folder,'mr_package_ivw.txt',sep='/')
  # unlink(fe_log_folder,recursive = TRUE)
  # con_fe <- file(fe_log_folder) # 
  # sink(con_fe, append=TRUE) # 
  # sink(con_fe, append=TRUE, type="message") # 记录message
  # print(model_result)
  # sink()
  # sink(type="message")
  # cat(readLines(fe_log_folder), sep="\n")
  # mr_re_I2<-read.delim(fe_log_folder,sep = '\t')
  # unlink(fe_log_folder,recursive = TRUE)
  # mr_re_I2<-mr_re_I2$Inverse.variance.weighted.method[nrow(mr_re_I2)]
  # mr_re_I2<-gsub("[[:space:]]", "", strsplit(mr_re_I2,split='=')[[1]][4])
  # mr_re_I2<-as.numeric(substr(mr_re_I2,1,nchar(mr_re_I2)-2))/100
  # if (is.na(mr_re_I2)){
  #   mr_re_I2<-'NA'
  # }  
  # all_res_forest1[all_res_forest1$Method=="MR IVW(RE)",'maf_I2']<-mr_re_I2
  # all_res_forest1[all_res_forest1$Method=="MR IVW(RE)",'r2_eaf']<-mr_re_I2
  # }

  
  F_r_stats=list('r2_maf_sum'=r2_maf_sum,'F_maf'=Fstat_maf,'Corhan_Q_I2_by_maf_result'=round((Q_ivw-Q_ivw_df)/Q_ivw,4),'r2_eaf_sum'=r2,'F_eaf'=Fstat)
  save(F_r_stats,file='F_r_stats.RData')
  
  write.csv(all_res_forest1,'original_har_radial_outliners_single_mr_result(sensitivity)_forest.csv', row.names = FALSE)
  write.csv(all_res,'original_har_radial_outliners_single_mr_result(sensitivity).csv', row.names = FALSE)
  whole_re<-rbind(whole_re,all_res)
  whole_re_forest<-rbind(whole_re_forest,all_res_forest1)
  head(all_res[,c("Method","outcome","exposure","nsnp","b","se","or", "or_lci95",  
                  "or_uci95","pval","intercept","intercept_se","intercept_pval",
                  "Q","Q_df","Q_pval","presso_pval")])
}


save(whole_re,file=paste(whole_result_folder,'original_har_radial_outliners_whole_mr_result(sensitivity).RData',sep='/'))
save(whole_re_forest,file=paste(whole_result_folder,'original_har_radial_outliners_whole_mr_result(sensitivity)_forest.RData',sep='/'))
write.csv(whole_re,paste(whole_result_folder,'original_har_radial_outliners_whole_mr_result(sensitivity).csv',sep='/'), row.names = FALSE)
write.csv(whole_re_forest,paste(whole_result_folder,'original_har_radial_outliners_whole_mr_result(sensitivity)_forest.csv',sep='/'), row.names = FALSE)
print(paste0('time spent:',' ',proc.time()-ptm))
table(whole_re$outcome,useNA = 'always')
table(whole_re$exposure,useNA = 'always')



