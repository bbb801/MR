rm(list=ls())
gc()
library("readxl")
library(TwoSampleMR)
library(tidyverse)
library(RadialMR)


server<-'home'
whole_result_folder<-'.../2_20/result'
cvs_path<-'.../data/data_table.csv'


d<-read_csv(cvs_path)
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

#write.csv(rbind(d_out,d_mcv,d_rdw),"data_table.csv", row.names = FALSE)

###################################################################################################
for (i in 1:nrow(d_out)){
  for (j in c('mcv','rcdw')){
    setwd(whole_result_folder)
    exp_name<-j
    check_fun(out_index=i,exp_name=exp_name)
  }
}
 i<-6
exp_name<-'mcv'
re<-check_fun(out_index=i,exp_name=exp_name)


check_fun<-function(out_index,exp_name){
  if((exp_name=='rcdw')|(exp_name=='rdw')){
    
    exp_id<-d_rdw$id#'
    exp_trait<-d_rdw$trait
  } else if(exp_name=='mcv'){
    exp_id<-d_mcv$id#
    exp_trait<-d_mcv$trait
  } else(stop('Check exposure!'))
  out_id<-d_out$id[i]#
  out_trait<-d_out$trait[i]
  save_folder<-paste(whole_result_folder,exp_id,out_id,sep='/')
  unlink(save_folder,recursive = TRUE)
  dir.create(save_folder, showWarnings = TRUE, recursive = TRUE)

  log_file<-paste(save_folder,paste(exp_trait,'_',out_trait,'pre.log',sep=''),sep='/')
  con <- file(log_file)
  sink(con, append=TRUE) # 
  sink(con, append=TRUE, type="message") # 
  print(paste('Exposure: ',exp_trait,'; ID ',exp_id))
  print(paste('Outcome: ',out_trait,'; ID ',out_id))
  
  if(!file.exists(paste(whole_result_folder,paste0(exp_id,'.RData'),sep='/'))){
    exp_original<-extract_instruments(
      outcomes=exp_id,
      clump = F
    )
    save(exp_original,file = paste(whole_result_folder,paste0(exp_id,'.RData'),sep='/'))
  }

  exp <- extract_instruments(
    outcomes=exp_id,                  
    clump=TRUE, #r2=0.01,
    kb=10000,access_token= NULL
  )  
  dput(names(exp))
  save(mcv1,file='...vegf.RData')
  
  
  out <- extract_outcome_data(
    snps=exp$SNP,
    outcomes=out_id,
    proxies = TRUE,
   # maf_threshold = 0.01,
    access_token = NULL
  )
  
  har <- harmonise_data(
    exposure_dat=exp,
    outcome_dat=out,
    action= 3)  
  

  res <- mr(har)
  res  #
  
  Mydata<-har
  Mydata<-Mydata[Mydata["ambiguous" ]==FALSE,]
  print(paste0('Ori SNP No.after deleting ambiguous:',nrow(Mydata)))
  
  dat<-Mydata
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
  print(paste0('Before Radial has ',nrow(Mydata),' snps; Then has ', nrow(dat2), ' snps.'))
  dat2<-steiger_filtering(dat2)
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
  print('*********************************************************************************')
  print(paste('Exposure: ',exp_trait,'; ID ',exp_id))
  print(paste('Outcome: ',out_trait,'; ID ',out_id))
  print('*********************************************************************************')
  print('Final Radial IVW*************************************************************************')
  ivwrad3 <- ivw_radial(raddat,0.05/nrow(raddat),3,summary=TRUE)
  print('Final Radial egger*************************************************************************')
  egger.model3<-egger_radial(raddat,0.05/nrow(raddat),3,summary=TRUE)
  het<-mr_heterogeneity(dat2, method_list=c("mr_egger_regression", "mr_ivw"))  
  print(rep('#',10))
  print('Heterogeneity test')  
  print(het)
  # 
  plt <- mr_pleiotropy_test(dat2) ####egger
  print(rep('#',10))
  print('Pleiotropy test')  
  print(plt)
  sink()
  sink(type="message")
  cat(readLines(log_file), sep="\n")
#  return(list_data <- list('exp_id'=exp_id,'exp_trait'=exp_trait,'out_id'=out_id,'out_trait'=out_trait,'ivw_result'=ivw_radial(raddat,0.05/nrow(raddat),3,summary=TRUE),'egger_result'=egger_radial(raddat,0.05/nrow(raddat),3,summary=TRUE)))
}  






