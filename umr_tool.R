
library("readxl")
library(TwoSampleMR)
library(phenoscanner)
library(LDlinkR)
library(dplyr)
library(cause)
library(RadialMR)
library(tidyverse)
library(data.table)
library(MRPRESSO)
library(MendelianRandomization)

d_fun<-function(data_list,data_table_csv,id=NA){
  data_list<-data.frame(data_list)
  data_list$folder<-dirname(data_list$data_list)
  
  
  data_list$id<-basename(unlist(lapply(strsplit(data_list$data_list,split='.',fixed=TRUE),head,1)))
  data_list$id<-gsub("finngen_R5_","finn-b-",data_list$id)
  data_list<-as_tibble(data_list)
  data_table<-read_csv(data_table_csv)
  final<-inner_join(data_list, data_table,#[,c('id','trait')], 
                    c("id" = "id"))
  if (!is.na(id)){
    final<-final[final$id%in%id,]
  }
  return(final)
}
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

filter_fun<-function(dat,folder,out_id,exp_id,out_trait,exp_trait,bfile,presso_needed=FALSE,save_ld=TRUE,pattern='merge_clump'){
  setwd(folder)
  if (pattern=='merge_clump'){
    raddat <- format_radial(BXG=dat$beta_hat_1, BYG=dat$beta_hat_2, 
                            seBXG=dat$seb1, seBYG=dat$seb2, RSID=dat$rsid)
  } else{
    dat<-steiger_filtering(dat)%>%subset(steiger_dir==TRUE)
    raddat <- format_radial(BXG=dat$beta.exposure, BYG=dat$beta.outcome, 
                            seBXG=dat$se.exposure, seBYG=dat$se.outcome, RSID=dat$SNP)
  }

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
      print(rep('2--',30))
    }
  }
  
  
  
  if (pattern=='merge_clump'){
    raddat$id.exposure<-exp_id
    raddat$id.outcome<-out_id
    raddat$exposure<-exp_trait
    raddat$outcome<-out_trait
    raddat$mr_keep<-TRUE
  } else {
    raddat<-dat%>%filter(SNP %in% raddat$SNP)
  }

  
  # print(paste0('Before Radial has ',nrow(dat),' snps; Then has ', nrow(raddat), ' snps.'))
  
  
  # print(paste0('After Steiger has ', nrow(dat2), ' snps.'))
  # 
  
  print(rep(0,20))
  mr_package_input_snp_corr<-ieugwasr::ld_matrix(
    raddat$SNP,
    plink_bin = genetics.binaRies::get_plink_binary(),
    bfile = bfile
  )
  print(rep(1,20))
  #mr_package_input_snp_corr<-ld_matrix(raddat$SNP, with_alleles = TRUE, pop = "EUR")
  if (save_ld==TRUE){
    save(mr_package_input_snp_corr,file='mr_package_input_snp_correlation.RData')
  }
  
  mr_package_input<-mr_input(bx = raddat$beta.exposure, bxse = raddat$se.exposure, by = raddat$beta.outcome, byse = raddat$se.outcome,snps = raddat$SNP,correlation=mr_package_input_snp_corr)
  
  # # rm(mr_package_input_snp_corr)
  #  gc()
  mr_package_lasso<-mr_lasso(mr_package_input)
  if (mr_package_lasso@SNPs!=mr_package_lasso@Valid){
    print(paste('Outcome: ',out_trait,'; ID ',out_id))
    stop('Lasso SNPS have outliers!')
  }
  
  
  singlesnp_results=mr_singlesnp(raddat, parameters = default_parameters(),
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
    raddat<-raddat[which(!(raddat$SNP %in% out_snp)),]
    
  } else if (all(single_snp_all$p<0.05)&any(single_snp_single$p>=0.05)){
    print('Outliners in single snp mr result:')
    print(singlesnp_results[single_snp_single$p>=0.05,][,c('SNP','p')])
    out_snp<-singlesnp_results[single_snp_single$p>=0.05,]$SNP
    raddat<-raddat[which(!(raddat$SNP %in% out_snp)),]
  }
  print(paste0('After single sno test has ', nrow(raddat), ' snps.'))
  rownames(raddat)<-NULL
  
  
  print('Final Radial IVW')
  ivwrad3 <- ivw_radial(raddat,0.05/nrow(raddat),3)
  print('Final Radial egger')
  egger.model3<-egger_radial(raddat,0.05/nrow(raddat),3)
  
  
  mr_results <- mr(raddat)
  res<-mr_results
  #l[[i]]<-res
  print(rep('#',10))
  print('MR result')
  print(res)
  
  het<-mr_heterogeneity(raddat, method_list=c("mr_egger_regression", "mr_ivw"))  #异质性检验——IVWorMRegger
  print(rep('#',10))
  print('Heterogeneity test')  
  print(het)
  # 
  plt <- mr_pleiotropy_test(raddat) ####egger
  print(rep('#',10))
  print('Pleiotropy test')  
  print(plt)
  
  if (presso_needed){
    presso<-mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data =raddat, NbDistribution = 5000, SignifThreshold = 0.05)
    print(rep('#',10))
    print('presso') 
    print(presso)
    save(presso,file='presso.RData')
  }
  
  return(raddat)
}
d_table_fun<-function(csv_path){
  d<-read_csv(csv_path)
  d_out<-d[grep(pattern="Aneurysms, operations, SAH|Subarachnoid haemmorrhage|Traumatic subarachnoid haemorrhage|Intracerebral haemmorrhage|Intracranial trauma|Nontraumatic intracranial haemmorrhage",d$trait),]
  d_out<-d_out[!grepl(pattern="no controls excluded|ICD10|Benign intracranial hypertension|Intracranial volume",d_out$trait),]
  d_out<-d_out[grep(pattern="European",d_out$population),]
  d_out<-d_out[!grepl(pattern="ukb",d_out$id),]
  d_out$trait
  d_out
  
  d_mcv<-d[grep(pattern="Mean corpuscular volume",d$trait),]
  d_mcv<-d_mcv[grep(pattern="European",d_mcv$population),]
  d_mcv<-d_mcv[grep(pattern="ukb",d_mcv$id),]
  d_mcv$trait
  d_mcv
  
  
  d_rdw<-d[grep(pattern="Red cell distribution width",d$trait),]
  d_rdw<-d_rdw[grep(pattern="European",d_rdw$population),]
  d_rdw$trait
  d_rdw
  return(list('out'=d_out,'mcv'=d_mcv,'rcdw'=d_rdw))
}


nn<-function (res, het, plt, sin, ao_slc = T, Exp = F, split.exposure = F, 
              split.outcome = F) 
  
{
  requireNamespace("plyr", quietly = TRUE)
  het <- het[, c("id.exposure", "id.outcome", "method", 
                 "Q", "Q_df", "Q_pval")]
  # Pos <- which(unlist(lapply(names(res), FUN = function(x) class(res[, 
  #                                                                    x]))) == "factor")
  # for (i in 1:length(Pos)) {
  #   res[, Pos[i]] <- as.character(res[, Pos[i]])
  # }
  # Pos <- which(unlist(lapply(names(het), FUN = function(x) class(het[, 
  #                                                                    x]))) == "factor")
  # for (i in 1:length(Pos)) {
  #   het[, Pos[i]] <- as.character(het[, Pos[i]])
  # }
  # Pos <- which(unlist(lapply(names(sin), FUN = function(x) class(sin[, 
  #                                                                    x]))) == "factor")
  # for (i in 1:length(Pos)) {
  #   sin[, Pos[i]] <- as.character(sin[, Pos[i]])
  # }
  sin <- sin[grep("[:0-9:]", sin$SNP), ]
  sin$method <- "Wald ratio"
  names(sin)[names(sin) == "p"] <- "pval"
  names(res)[names(res) == "method"] <- "Method"
  names(het)[names(het) == "method"] <- "Method"
  names(sin)[names(sin) == "method"] <- "Method"
  res <- merge(res, het, by = c("id.outcome", "id.exposure", 
                                "Method"), all.x = T)
  res <- plyr::rbind.fill(res, sin[, c("exposure", "outcome", 
                                       "id.exposure", "id.outcome", "SNP", 
                                       "b", "se", "pval", "Method")])
  if (ao_slc) {
    ao <- available_outcomes(access_token = ieugwasr::check_access_token())
    names(ao)[names(ao) == "nsnp"] <- "nsnps.outcome.array"
    res <- merge(res, ao[, !names(ao) %in% c("unit", 
                                             "priority", "sd", "path", "note", 
                                             "filename", "access", "mr")], by.x = "id.outcome", 
                 by.y = "id")
  }
  res$nsnp[is.na(res$nsnp)] <- 1
  for (i in unique(res$id.outcome)) {
    Methods <- unique(res$Method[res$id.outcome == i])
    Methods <- Methods[Methods != "Wald ratio"]
    for (j in unique(Methods)) {
      res$SNP[res$id.outcome == i & res$Method == j] <- paste(res$SNP[res$id.outcome == 
                                                                        i & res$Method == "Wald ratio"], collapse = "; ")
    }
  }
  if (Exp) {
    res$or <- exp(res$b)
    res$or_lci95 <- exp(res$b - res$se * 1.96)
    res$or_uci95 <- exp(res$b + res$se * 1.96)
  }
  plt <- plt[, c("id.outcome", "id.exposure", "egger_intercept", 
                 "se", "pval")]
  plt$Method <- "MR Egger"
  names(plt)[names(plt) == "egger_intercept"] <- "intercept"
  names(plt)[names(plt) == "se"] <- "intercept_se"
  names(plt)[names(plt) == "pval"] <- "intercept_pval"
  res <- merge(res, plt, by = c("id.outcome", "id.exposure", 
                                "Method"), all.x = T)
  if (split.exposure) {
    res <- split_exposure(res)
  }
  if (split.outcome) {
    res <- split_outcome(res)
  }
  Cols <- c("Method", "outcome", "exposure", 
            "nsnp", "b", "se", "pval", "intercept", 
            "intercept_se", "intercept_pval", "Q", 
            "Q_df", "Q_pval", "consortium", "ncase", 
            "ncontrol", "pmid", "population")
  res <- res[, c(names(res)[names(res) %in% Cols], names(res)[which(!names(res) %in% 
                                                                      Cols)])]
  return(res)
}
p2n<-function(model_result,folder){
  fe_log_folder<-paste(folder,'mr_package_ivw.txt',sep='/')
  unlink(fe_log_folder,recursive = TRUE)
  con_fe <- file(fe_log_folder) # 
  sink(con_fe, append=TRUE) # 
  sink(con_fe, append=TRUE, type="message") # 
  print(model_result)
  sink()
  sink(type="message")
  cat(readLines(fe_log_folder), sep="\n")
  a<-read.delim(fe_log_folder,sep = '\t')
  unlink(fe_log_folder,recursive = TRUE)
  a<-a$Inverse.variance.weighted.method[nrow(a)]
  a<-gsub("[[:space:]]", "", strsplit(a,split='=')[[1]][4])
  a<-as.numeric(substr(a,1,nchar(a)-2))/100
  if (is.na(a)){
    return('NA')
  } else {return(a)}
}

p2n2<-function(model_result,folder){
  fe_log_folder<-paste(folder,'mr_package_ivw.Rout',sep='/')
  unlink(fe_log_folder,recursive = TRUE)
  con_fe <- file(fe_log_folder,open = 'wt') # 
  sink(con_fe) # 
  sink(con_fe, type="message") # 
  print(model_result)
  sink(type="message")
  sink()
  file.show(fe_log_folder)
  a<-read.delim(fe_log_folder,sep = '\t')
  unlink(fe_log_folder,recursive = TRUE)
  a<-a$Inverse.variance.weighted.method[nrow(a)]
  a<-gsub("[[:space:]]", "", strsplit(a,split='=')[[1]][4])
  a<-as.numeric(substr(a,1,nchar(a)-2))/100
  return(a)
}



or_fun_two_step<-function (data1,all_res_forest,folder,raddat,type='presso') 
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
    con_ivw <- file(ivw_log_folder) # 
    sink(con_ivw, append=TRUE) # 
    sink(con_ivw, append=TRUE, type="message") # 
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
    F_value<-as.numeric(str_extract_all(ivw_txt$Radial.IVW[7],"\\d+\\.\\d+|\\d+")[[1]])[1]
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
    mr_res$F_value<-F_value
    mr_res[,c("intercept", "intercept_se", "intercept_pval",'SNP','presso_pval')]<-NA
  }
  else  if (type=='rad_egger'){
    egger_log_folder<-paste(folder,'rad_egger.txt',sep='/')
    unlink(egger_log_folder,recursive = TRUE)
    con_egger <- file(egger_log_folder) # 
    sink(con_egger, append=TRUE) # 
    sink(con_egger, append=TRUE, type="message") # 
    egger_radial(raddat,0.05/nrow(raddat),3)
    sink()
    sink(type="message")
    cat(readLines(egger_log_folder), sep="\n")
    egger_txt<-read.delim(egger_log_folder,sep = "\t")
    unlink(egger_log_folder,recursive = TRUE)
    egger_Q<-as.numeric(str_extract_all(egger_txt$Radial.MR.Egger[6],"\\d+\\.\\d+|\\d+")[[1]])[1]
    egger_Q_df<-as.numeric(str_extract_all(egger_txt$Radial.MR.Egger[6],"\\d+\\.\\d+|\\d+")[[1]])[2]
    egger_Q_p<-as.numeric(str_extract_all(egger_txt$Radial.MR.Egger[6],"\\d+\\.\\d+|\\d+")[[1]])[3]
    F_value<-as.numeric(str_extract_all(egger_txt$Radial.MR.Egger[5],"\\d+\\.\\d+|\\d+")[[1]])[1]
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
    mr_res$F_value<-F_value
    
    mr_res[,c('SNP')]<-NA
    mr_res$presso_pval<-NA
  } else if (type=='mr_package_ivw'){
    mr_res<-all_res_forest[1:2,]
    mr_res$Method<-c('MR IVW(FE)','MR IVW(RE)')
    mr_res$nsnp<-c(data1@SNPs,mr_package_ivw_re@SNPs)
    mr_res$b<-c(data1@Estimate,mr_package_ivw_re@Estimate)
    
    mr_res$se<-c(data1@StdError,mr_package_ivw_re@StdError)
    mr_res$pval<-c(data1@Pvalue,mr_package_ivw_re@Pvalue)
    mr_res$lo_ci <-c(data1@CILower,mr_package_ivw_re@CILower)
    mr_res$up_ci <-c(data1@CIUpper,mr_package_ivw_re@CIUpper)
    mr_res$or <- exp(mr_res$b)
    mr_res$or_lci95 <- exp(mr_res$lo_ci)
    mr_res$or_uci95 <- exp(mr_res$up_ci)
    mr_res$lo_ci<-NULL
    mr_res$up_ci<-NULL
    mr_res$nsnp<-c(data1@SNPs,mr_package_ivw_re@SNPs)
    mr_res$Q<-c(data1@Heter.Stat[1],mr_package_ivw_re@Heter.Stat[1])
    mr_res$Q_pval<-c(data1@Heter.Stat[2],mr_package_ivw_re@Heter.Stat[2])
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

#    }