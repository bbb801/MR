library(TwoSampleMR)
library("readxl")
library(tidyverse)
library(RadialMR)
library(boot)
library(LDlinkR)
library(dplyr)
library(MVMR)
library(MRPRESSO)
library(boot.pval)
library(MendelianRandomization)
rm(list=ls())
gc()
server<-'home'
mvmr_tool<-'.../2_20/code/mvmr_tool.R'
whole_result_folder<-'.../two_step'
b1_DBP<-'.../Action_3_after_lasso_har_data.RData'#mcv RCDW & diastolic blood pressure
b1_SBP<-'.../ieu-b-38/Action_3_after_lasso_har_data.RData'#mcv RCDW & diastolic blood pressure
b0_ICH<-'...//finn-b-I9_ICH/Action_3_after_lasso_har_data.RData'
b2_DBP_ICH<-'.../RCDW_MCV_DBP/finn-b-I9_ICH/Action_3_after_lasso_har_data.RData'
b2_SBP_ICH<-'...//RCDW_MCV_SBP/finn-b-I9_ICH/Action_3_after_lasso_har_data.RData'
b0_point<-'.../RCDW_MCV/whole_result_action3_new.RData'
b1_point<-'.../whole_result_action3_new.RData'
b2_point<-'.../whole_result_action3_new.RData'
b0_nITH<-'.../finn-b-I9_INTRACRA/Action_3_after_lasso_har_data.RData'
b2_DBP_nITH<-'.../finn-b-I9_INTRACRA/Action_3_after_lasso_har_data.RData'
b1_DBP_exp<-'.../exposure_dat3.RData'


source(mvmr_tool)


###setting


############################################

mvmr_fe_b <- function(data, indices,exp_abbr,out_abbr,stats,logv) { 
  d1 <- data[indices,]
  bx_matrix = as.matrix(as.matrix(dplyr::select(d1,ends_with('_b'),-paste0(out_abbr,'_b'))))
  sebx_matrix = as.matrix(as.matrix(dplyr::select(d1,ends_with('_se'),-paste0(out_abbr,'_se'))))
  
  mvmr_dat_input<-mr_mvinput(bx = bx_matrix, bxse = sebx_matrix,
                             by = d1[,paste0(out_abbr,'_b')],
                             byse = d1[,paste0(out_abbr,'_se')],snps =d1$SNP,outcome = out_name_new_list(out_abbr,id2abbr_real=1,reverse=TRUE,style='American'),exposure = out_name_new_list(exp_abbr,id2abbr_real=1,reverse=TRUE,style='American'),effect_allele =d1$allele_0,other_allele =d1$allele_1  )
  mvmr_ivw_fe = mr_mvivw(mvmr_dat_input,model='fixed')
  
  if (stats=='b2'){
    rn<-rownames(as.data.frame(mvmr_ivw_fe@Estimate))
    med_b<-setdiff(rn,c("BxRCDW_b" ,"BxMCV_b" ))
    bb<-as.data.frame(mvmr_ivw_fe@Estimate)[med_b,]
  } else{
    bb<-as.data.frame(mvmr_ivw_fe@Estimate)['BxMCV_b',]
  }
  if (logv){
    return(exp(bb))
  } else{return(bb)}
  
}

boot_fun<-function(stats='b0',med='DBP',out='ICH',boot_n=boot_n,logv=TRUE){
  ptm <- proc.time()
  if(stats=='b0'){
    folder<-paste(whole_result_folder,stats,out,sep='/')
    if (out=='ICH'){
      load(b0_ICH)
      bd<-rawdat_mvmr
      rm(rawdat_mvmr)
      load(b0_point)
      point_dat<-whole_result%>%filter(Method == "MR IVW(FE)",id.outcome=="finn-b-I9_ICH",id.exposure=='ukb-d-30040_irnt')
      out_abbr<-out
      exp_abbr<-head(gsub('_se','',names(bd)[endsWith(names(bd), '_se')]),-1)
      
    }else if (out=='nITH'){
      load(b0_nITH)
      bd<-rawdat_mvmr
      rm(rawdat_mvmr)
      load(b0_point)
      point_dat<-whole_result%>%filter(Method == "MR IVW(FE)",id.outcome=="finn-b-I9_INTRACRA",id.exposure=='ukb-d-30040_irnt')
      out_abbr<-out
      exp_abbr<-head(gsub('_se','',names(bd)[endsWith(names(bd), '_se')]),-1)
      
    }
  }else if(stats=='b1'){
    folder<-paste(whole_result_folder,stats,med,sep='/')
    load(b1_point)
    if (med=='SBP'){
      load(b1_SBP)
      bd<-rawdat_mvmr
      rm(rawdat_mvmr)
      point_dat<-whole_result%>%filter(Method == "MR IVW(FE)",id.outcome=="ieu-b-38",id.exposure=='ukb-d-30040_irnt')
    }else if(med=='DBP'){
      load(b1_DBP)
      bd<-rawdat_mvmr
      rm(rawdat_mvmr)
      point_dat<-whole_result%>%filter(Method == "MR IVW(FE)",id.outcome=="ieu-b-39",id.exposure=='ukb-d-30040_irnt')
    } 
    out_abbr<-med
    exp_abbr<-head(gsub('_se','',names(bd)[endsWith(names(bd), '_se')]),-1)
  }else if((stats=='b2')|(stats=='b0_direct')){
    folder<-paste(whole_result_folder,stats,paste0('MCV_',med),out,sep='/')
    load(b2_point)
    if (med=='SBP'){
      if(out=='ICH'){
        load(b2_SBP_ICH)
        bd<-rawdat_mvmr
        rm(rawdat_mvmr)
        point_dat_mcv<-whole_result%>%filter(Method == "MR IVW(FE)",id.outcome=="finn-b-I9_ICH",id.exposure=='ukb-d-30040_irnt')%>%filter(str_detect(exp_group,med))
        point_dat_med<-whole_result%>%filter(Method == "MR IVW(FE)",id.outcome=="finn-b-I9_ICH",id.exposure=='ieu-b-38')%>%filter(str_detect(exp_group,med))
      }else if(out=='nITH'){
        load(b2_SBP_nITH)
        bd<-rawdat_mvmr
        rm(rawdat_mvmr)
        point_dat_mcv<-whole_result%>%filter(Method == "MR IVW(FE)",id.outcome=="finn-b-I9_INTRACRA",id.exposure=='ukb-d-30040_irnt')%>%filter(str_detect(exp_group,med))
        point_dat_med<-whole_result%>%filter(Method == "MR IVW(FE)",id.outcome=="finn-b-I9_INTRACRA",id.exposure=='ieu-b-38')%>%filter(str_detect(exp_group,med))
      }
      
    }else if(med=='DBP'){
      if(out=='ICH'){
        load(b2_DBP_ICH)
        bd<-rawdat_mvmr
        rm(rawdat_mvmr)
        #b2_point
        point_dat_mcv<-whole_result%>%filter(Method == "MR IVW(FE)",id.outcome=="finn-b-I9_ICH",id.exposure=='ukb-d-30040_irnt')%>%filter(str_detect(exp_group,med))
        point_dat_med<-whole_result%>%filter(Method == "MR IVW(FE)",id.outcome=="finn-b-I9_ICH",id.exposure=='ieu-b-39')%>%filter(str_detect(exp_group,med))
      }else if(out=='nITH'){
        load(b2_DBP_nITH)
        bd<-rawdat_mvmr
        rm(rawdat_mvmr)
        point_dat_mcv<-whole_result%>%filter(Method == "MR IVW(FE)",id.outcome=="finn-b-I9_INTRACRA",id.exposure=='ukb-d-30040_irnt')%>%filter(str_detect(exp_group,med))
        point_dat_med<-whole_result%>%filter(Method == "MR IVW(FE)",id.outcome=="finn-b-I9_INTRACRA",id.exposure=='ieu-b-39')%>%filter(str_detect(exp_group,med))
      }

    }
    out_abbr<-out
    exp_abbr<-head(gsub('_se','',names(bd)[endsWith(names(bd), '_se')]),-1)
    if(stats=='b2'){
      point_dat<-point_dat_med
    }else if(stats=='b0_direct'){
      point_dat<-point_dat_mcv
    }
  }
  
  if(logv==TRUE){
    point_dat<-point_dat%>%subset(select=c('or_lci95','or',  'or_uci95','pval'))
  } else {
    point_dat<-point_dat%>%subset(select=c('lo_ci','b',  'up_ci','pval'))
  }

  gc()

  
  
  dir.create(folder, showWarnings = TRUE, recursive = TRUE)
  setwd(folder)
  set.seed(1234)
  boot_b<-boot(bd, mvmr_fe_b,R=boot_n,exp_abbr=exp_abbr,out_abbr=out_abbr,stats=stats,logv=logv)

  save(boot_b,file=paste(folder,paste0('b_',boot_n,'.RData' ),sep='/'))
  res<-list('boot_b'=boot_b,'point_dat'=point_dat)
   # boot_list<-list('b0'=boot_b0,'b1'=boot_b1,'b2'=boot_b2)

  print(paste0('time spent:',' ',proc.time()-ptm))
  print('Done!!!!')
  return(res)

}
#DBP
#diff
develop<-FALSE
model<-'FE'#MVMR
med<-'DBP'
boot_n<-1000
out<-'ICH'
b0<-boot_fun(stats='b0',out=out,boot_n=boot_n,logv=FALSE)
b1<-boot_fun(stats='b1',med=med,boot_n=boot_n,logv=FALSE)
b2<-boot_fun(stats='b2',med=med,out=out,boot_n=boot_n,logv=FALSE)
b0_direct<-boot_fun(stats='b0_indirect',med=med,out=out,boot_n=boot_n,logv=FALSE)

b0_exp<-boot_fun(stats='b0',out=out,boot_n=boot_n,logv=TRUE)
# 
# t_indirect_diff<-b0$boot_b$t-b0_direct$boot_b$t
# t.test(t_indirect_diff, mu = 0, alternative = "two.sided")
# 
# proportion_diff<-(b0$boot_b$t-b0_direct$boot_b$t)/b0$boot_b$t
# quantile(proportion_diff,c(0.25,0.75))
# proportion_diff_point<-(b0$point_dat$b-b0_direct$point_dat$b)/b0$point_dat$b
# ####
# #product
# t_indirect_product<-b1$boot_b$t*b2$boot_b$t
# t.test(t_indirect_product, mu = 0, alternative = "two.sided")
# 
# proportion_product<-(t_indirect_product)/b0$boot_b$t
# quantile(proportion_product,c(0.25,0.75))
# proportion_product_point<-(b1$point_dat$b*b2$point_dat$b)/b0$point_dat$b

#proportion-mediated measure by point

b0_direct_exp<-boot_fun(stats='b0_direct',med=med,out=out,boot_n=boot_n,logv=TRUE)
#t.test(b0_direct_exp$boot_b$t, mu = 1, alternative = "two.sided")
b0_direct_exp$point_dat$or*(exp(b1$point_dat$b*b2$point_dat$b)-1)/(b0_direct_exp$point_dat$or*exp(b1$point_dat$b*b2$point_dat$b)-1)# proportion-mediated measure by point
pm<-b0_direct_exp$boot_b$t*(exp(b1$boot_b$t*b2$boot_b$t)-1)/(b0_direct_exp$boot_b$t*exp(b1$boot_b$t*b2$boot_b$t)-1)
quantile(pm,c(0.025,0.975))

#indirect
exp(b1$point_dat$b*b2$point_dat$b)
t_indirect_product_exp<-exp(b1$boot_b$t*b2$boot_b$t)
t.test(t_indirect_product_exp, mu = 1, alternative = "two.sided")
quantile(t_indirect_product_exp,c(0.025,0.975))



product_fun<-function(med='DBP',out='ICH'){
  b0<-boot_fun(stats='b0',out=out,boot_n=boot_n,logv=FALSE)
  b1<-boot_fun(stats='b1',med=med,boot_n=boot_n,logv=FALSE)
  b2<-boot_fun(stats='b2',med=med,out=out,boot_n=boot_n,logv=FALSE)
  b0_direct<-boot_fun(stats='b0_direct',med=med,out=out,boot_n=boot_n,logv=FALSE)
  b0_exp<-boot_fun(stats='b0',out=out,boot_n=boot_n,logv=TRUE)
  b0_direct_exp<-boot_fun(stats='b0_direct',med='DBP',out='ICH',boot_n=boot_n,logv=TRUE)
  print('product of coefficients method: Exp stats')
  product_indirect_est<-exp(b1$point_dat$b*b2$point_dat$b)
  t_indirect_product_exp<-exp(b1$boot_b$t*b2$boot_b$t)
  product_indirect_test<-t.test(t_indirect_product_exp, mu = 1, alternative = "two.sided")
  product_indirect_ci<-quantile(t_indirect_product_exp,c(0.025,0.975))
  print(paste0('Indirect effect estimates on an odds ratio scale: ',product_indirect_est))
  print('t-test for indirect effect estimates on an odds ratio scale: ')
  print(product_indirect_test)
  print('95% CI for indirect effect estimates on an odds ratio scale: ')
  print(product_indirect_ci)
  
  pm_est<-b0_direct_exp$point_dat$or*(exp(b1$point_dat$b*b2$point_dat$b)-1)/(b0_direct_exp$point_dat$or*exp(b1$point_dat$b*b2$point_dat$b)-1)# proportion-mediated measure by point
  pm<-b0_direct_exp$boot_b$t*(exp(b1$boot_b$t*b2$boot_b$t)-1)/(b0_direct_exp$boot_b$t*exp(b1$boot_b$t*b2$boot_b$t)-1)
  pm_test<-t.test(pm, mu = 0, alternative = "greater")
  pm_ci<-quantile(pm,c(0.025,0.975))
  print(paste0('Proportion-mediated measure: ',pm_est))
  print('95% CI for Proportion-mediated measure:')
  print(pm_ci)
  
  fram<-c('Exposure',	'Mediator',	'Outcome',	'Total_effect',	'Indirect_effect',	'Direct_effect',	'Proportion_mediated')
  rr<-data.frame(matrix(1:2*length(fram), nrow=1, ncol=length(fram)))
  names(rr)<-fram
  rr<-rr%>%mutate(EXposure='MCV',Mediator=med,Outcome=out,	Total_effect=paste0(round(b0_exp$point_dat$or,3),' (',round(b0_exp$point_dat$or_lci95,3),',',round(b0_exp$point_dat$or_uci95,3),')', '; P=',round(b0_exp$point_dat$pval,3)),
                  Indirect_effect=paste0(round(product_indirect_est,3),' (',round(product_indirect_ci[1],3),',',round(product_indirect_ci[2],3),')', '; P=',round(product_indirect_test$p.value,3)),
                  Direct_effect=paste0(round(b0_direct_exp$point_dat$or,3),' (',round(b0_direct_exp$point_dat$or_lci95,3),',',round(b0_direct_exp$point_dat$or_uci95,3),')', '; P=',round(b0_direct_exp$point_dat$pval,3)),
                  Proportion_mediated=paste0(round(pm_est,3),' (',round(pm_ci[1],3),round(pm_ci[2],3),')', '; P=',round(pm_test$p.value,3)))
  return(rr)

}

mediation_result<-rbind(product_fun(med='DBP',out='ICH'),product_fun(med='DBP',out='nITH'))
write.csv(mediation_result,".../hypertension/two_step/result.csv",row.names = F)



#cal r & F for overlap bias
n=mean(dat$samplesize.exposure)
k=13
p=dat$eaf.exposure
r2.1<-(2*(dat$beta.exposure^2)*p*(1-p))/1
r2=sum(r2.1)
r2 
Fstat <- r2*(n-1-k)/((1-r2)*k)
Fstat

load(b1_DBP)
final_har<-rawdat_mvmr
rm(rawdat_mvmr)
b1_DBP_exp<-'.../hypertension/systolic/mvmr_no_corr/RCDW_MCV/ieu-b-39/exposure_dat3.RData'
load(b1_DBP_exp)
exp_mvmr<-exposure_dat
rm(exposure_dat)

har_mvmr<-exp_mvmr%>%filter(id.exposure=='ukb-d-30040_irnt',SNP%in%final_har$SNP)

exp_umr <- extract_instruments(
  outcomes='ukb-d-30040_irnt',                  
  clump=FALSE, r2=0.01,
  kb=10000,access_token= NULL
)

har_umr<-exp_umr%>%filter(id.exposure=='ukb-d-30040_irnt',SNP%in%final_har$SNP)

n=350473
k=nrow(har_mvmr)
p=har_mvmr$eaf.exposure
r2.1<-(2*(har_mvmr$beta.exposure^2)*p*(1-p))/1
r2=sum(r2.1)
r2
Fstat <- r2*(n-1-k)/((1-r2)*k)
Fstat
n_out<-757601
