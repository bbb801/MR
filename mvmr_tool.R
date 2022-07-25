
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
library(tidyverse)
library(RColorBrewer)
library(cowplot)
library(gridExtra)
library(MASS)
library(glmnet)
library(quantreg)
library(robustbase)
library(Hmisc)

change_name_mvmr<-function(nam,id_outcome){
  n<-length(nam)/2
  for (i in 1:length(nam)){
    
    if (i<=n){
      if (nam[i]==out_name(nam[i])){
        outcome_abbr<-out_name(id_outcome)
        nam[i]<-paste0(outcome_abbr,'_b')
      } else{
        nam[i]<-paste0(out_name(nam[i]),'_b')
      }
    } else {
      if (nam[i]==out_name(nam[i])){
        outcome_abbr<-out_name(id_outcome)
        nam[i]<-paste0(outcome_abbr,'_se')
      }else{
        nam[i]<-paste0(out_name(nam[i]),'_se')
      }
    }
  }
  re<-list('renam'=nam,'out_abbr'=outcome_abbr)
  return(re)
}
clean_name_fun<-function(x,style='American'){
  x<-do.call('rbind',strsplit(x,split = '||',fixed = TRUE))[,1]
  x<-trimws(x, which = c("both", "left", "right"), whitespace = "[ \t\r\n]")
  if (style=='American'){
    x<-gsub("haemmorrhage|haemorrhage", "hemorrhage", x) 
    x<-gsub("SAH", "subarachnoid hemorrhage", x)
    
  } else{
    x<-gsub("haemmorrhage", "haemorrhage", x) 
    x<-gsub("SAH", "subarachnoid haemorrhage", x)
  }
  return(capitalize(x))
}

rbind_fun<-function(x,n){
  emp<-data.frame()
  for (i in 1:n){
    emp<-rbind(emp,x)
  }
  rownames(emp)<-NULL
  return(emp)
}
out_name_list<-function(x,id2abbr=TRUE){
  emp<-data.frame()
  for (i in x){
    if (id2abbr){
      emp<-c(emp,out_name(i))
    } else{
      emp<-c(emp,out_name2(i))
    }
    
  }
  return(unlist(emp))
}
out_name_new_list<-function(x,id2abbr_real=1,reverse=FALSE,style='American'){
  emp<-data.frame()
  for (i in x){

      emp<-c(emp,out_name_new(i,id2abbr_real=id2abbr_real,reverse=reverse,style=style))
    }
  return(unlist(emp))
}
out_name_new<-function(x,id2abbr_real=1,reverse=FALSE,style='American'){
  if (id2abbr_real==1&reverse==FALSE){
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
    }else if(x=='finn-b-I9_SAH'){
      return('SAH')
    }else if(x=='finn-b-I9_SAHANEUR'){
      return('AOS')
    }else{
      return(x)
    }
  } else if(id2abbr_real==1&reverse==TRUE){
    if(x=='RCDW'){
      return('ebi-a-GCST006804')
    } else if (x=='MCV'){
      return('ukb-d-30040_irnt')
    } else if(x=='SBP')
    {
      return('ieu-b-38')
    }else if(x=='DBP')
    {
      return('ieu-b-39')
    }else if(x=='BD')
    {
      return('finn-b-R18_ELEVATED_ERYTHROCYTE_SEDIM_RATE_ABNORMALITY_PLASMA_VISCO')
    }else if(x=='ICH')
    {
      return('finn-b-I9_ICH')
    }else if(x=='nITH')
    {
      return('finn-b-I9_INTRACRA')
    }else if(x=='SAH')
    {
      return('finn-b-I9_SAH')
    }else if(x=='AOS')
    {
      return('finn-b-I9_SAHANEUR')
    } else{
      return(x)
    }
  } else if(id2abbr_real==2&reverse==FALSE){
    d<-read.csv('.../2_20/alldata.csv', encoding="UTF-8")
    x<-d[d$id==x,]$trait
    x<-clean_name_fun(x,style=style)
    return(x)
  } else if(id2abbr_real==2&reverse==TRUE){
    d<-read.csv('.../2_20/alldata.csv', encoding="UTF-8")
    x<-d[d$trait==x,]$id
    return(x)
  }
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

out_name2<-function(x){
  if(x=='RCDW'){
    return('ebi-a-GCST006804')
  } else if (x=='MCV'){
    return('ukb-d-30040_irnt')
  } else if
  (x=='SBP'){
    return('ieu-b-38')
  }else if
  (x=='DBP'){
    return('ieu-b-39')
  }else if
  (x=='BD'){
    return('finn-b-R18_ELEVATED_ERYTHROCYTE_SEDIM_RATE_ABNORMALITY_PLASMA_VISCO')
    
  }else if
  (x=='ICH'){
    return('finn-b-I9_ICH')
  }else if(x=='nITH'){
    return('finn-b-I9_INTRACRA')
  } else{
    return(x)
  }
}

qhet_mvmr2<-function(r_input,pcor,CI,iterations){
  
  # convert MRMVInput object to mvmr_format
  if ("MRMVInput" %in% class(r_input)) {
    r_input <- mrmvinput_to_mvmr_format(r_input)
  }
  
  # Perform check that r_input has been formatted using format_mvmr function
  if(!("mvmr_format" %in%
       class(r_input))) {
    stop('The class of the data object must be "mvmr_format", please resave the object with the output of format_mvmr().')
  }
  
  requireNamespace("boot", quietly = TRUE)
  
  warning("qhet_mvmr() is currently undergoing development.")
  
  if(missing(CI)) {
    CI<-F
    warning("95 percent confidence interval not calculated")
  }
  
  if(missing(iterations)) {
    iterations<-1000
    warning("Iterations for bootstrap not specified. Default = 1000")
  }
  
  exp.number<-length(names(r_input)[-c(1,2,3)])/2
  
  Qtemp<-function(r_input,pcor){
    
    exp.number<-length(names(r_input)[-c(1,2,3)])/2
    stderr<-as.matrix(r_input[,(exp.number + 4):length(r_input)])
    correlation<-pcor
    gammahat<-r_input$betaYG
    segamma<-r_input$sebetaYG
    pihat<-as.matrix(r_input[,c(4:(3+exp.number))])
    
    #estimation with the heterogeneity statistic
    
    PL_MVMR = function(a){
      tau2   = a[1]
      
      PL2_MVMR = function(ab){
        b<-ab
        
        cov = matrix(nrow = exp.number, ncol = exp.number)
        w=NULL
        for(l in 1:nrow(r_input)){
          for(pp in 1:exp.number){
            for(p2 in 1:exp.number){
              cov[pp,p2] <- correlation[pp,p2]*stderr[l,pp]*stderr[l,p2]
            }
          }
          
          segamma<-r_input$sebetaYG
          
          w[l] <- segamma[l]^2+t(b)%*%cov%*%b + tau2
        }
        
        q =  sum((1/w)*((gammahat - pihat%*%b)^2))
        
        return(q)
        
      }
      
      st_PL2  =  rep(0,exp.number)
      
      bc    = optim(st_PL2, PL2_MVMR)
      
      bcresults <- bc$par
      
      cov = matrix(nrow =exp.number, ncol = exp.number)
      w=NULL
      for(l in 1:nrow(r_input)){
        for(pp in 1:exp.number){
          for(p2 in 1:exp.number){
            cov[pp,p2] <- correlation[pp,p2]*stderr[l,pp]*stderr[l,p2]
          }
        }
        
        w[l] <- segamma[l]^2+t(bcresults)%*%cov%*%bcresults + tau2
      }
      
      q = (sum((1/w)*((gammahat - pihat%*%bcresults)^2))-(nrow(r_input)-2))^2
      
    }
    
    PL2_MVMR = function(ab){
      b=ab
      
      w=NULL
      cov = matrix(nrow = exp.number, ncol = exp.number)
      for(l in 1:nrow(r_input)){
        for(pp in 1:exp.number){
          for(p2 in 1:exp.number){
            cov[pp,p2] <- correlation[pp,p2]*stderr[l,pp]*stderr[l,p2]
          }
        }
        
        w[l] <- segamma[l]^2+t(b)%*%cov%*%b + tau_i
      }
      
      q =  sum((1/w)*((gammahat - pihat%*%b)^2))
      
    }
    
    limltauest  = optimize(PL_MVMR,interval=c(-10,10))
    tau_i = limltauest$objective
    
    tau = tau_i
    
    liml_het2<- optim(rep(0,exp.number), PL2_MVMR)
    limlhets <- liml_het2$par
    Qexact_het <- liml_het2$value
    
    Effects<-limlhets
    Effects<-data.frame(Effects)
    names(Effects)<-"Effect Estimates"
    for(i in 1:exp.number){
      rownames(Effects)[i]<-paste("Exposure",i,sep=" ")
      
    }
    
    return(Effects)
    
    
  }
  
  if(CI==F){
    
    res<-Qtemp(r_input,pcor)
    
  }
  
  if(CI==T){
    
    bootse<-function(data, indices){
      
      bres<-Qtemp(data[indices,],pcor)[,1]
      
      return(bres)
      
    }
    
    b.results<-boot::boot(data=r_input,statistic=bootse,R=iterations)
    
    lcb<-NULL
    ucb<-NULL
    ci<-NULL
    
    
    for(i in 1:exp.number){
      
      lcb[i]<-round(boot::boot.ci(b.results, type="bca",index=i)$bca[4],digits=3)
      ucb[i]<-round(boot::boot.ci(b.results, type="bca",index=i)$bca[5],digits=3)
      
      ci[i]<-paste(lcb[i],ucb[i],sep="-")
      
    }
    
    res<-data.frame(b.results$t0,ci)
    
    names(res)<-c("Effect Estimates","95% CI")
    for(i in 1:exp.number){
      rownames(res)[i]<-paste("Exposure",i,sep=" ")
      
    }
    
  }
  
  return(list(res=res, b_results=b.results,boot_ci=boot::boot.ci(b.results, type="bca")))
  
  
}



mvmr_med_boot = function(bx, sebx, by, seby, N){
  est = sapply(1:N, function(i){
    p = length(by)
    k = dim(bx)[2]
    Sx = lapply(1:p, function(j){diag(sebx[j, ]^2)})
    bxboot = sapply(1:p, function(j){mvrnorm(1, bx[j, ], Sx[[j]])})
    bxboot = t(bxboot)
    byboot = rnorm(p, by, seby)
    rq(byboot ~ bxboot - 1, weights = seby^-2)$coefficients
  })
  apply(est, 1, sd)
}

p_fun<-function(b,se){
  nb<-length(b)
  pval<-pnorm(abs(b/se),lower.tail = F)*2
  names(pval)<-paste0('px',1:nb)
  return(pval)
}

mvmr_median = function(bx, sebx, by, seby, boot = FALSE, boot_it = 1000){
  qr_mod = rq(by ~ bx - 1, weights = seby^-2)
  if (boot == TRUE){
    boot_se = mvmr_med_boot(bx, sebx, by, seby, boot_it)
    pval<-p_fun(qr_mod$coefficients,boot_se)
    return(list("coefficients" = qr_mod$coefficients, "se" = boot_se,'pval'=pval))
  } else {
    return(list("coefficients" = qr_mod$coefficients))
  }
}
#mvmr_median(bx=bx_matrix, sebx=sebx_matrix, by=beta_se_all$Stroke_b, seby=beta_se_all$Stroke_se,boot=TRUE)

cv.mvmr_lasso = function(bx, by, seby){
  p = dim(bx)[1]
  k = dim(bx)[2]
  S = diag(seby^-2)
  b = S^(1/2) %*% bx
  Pb = b %*% solve(t(b) %*% b, t(b))
  xlas = (diag(p) - Pb) %*% S^(1/2)
  ylas = (diag(p) - Pb) %*% S^(1/2) %*% by
  alas = glmnet(xlas, ylas, intercept = FALSE)
  lamseq = sort(alas$lambda)
  lamlen = length(lamseq)
  rse = sapply(1:lamlen, function(j){
    av = which(alas$beta[, (lamlen - j + 1)] == 0)
    mod = lm.fit(as.matrix(S[av, av]^(1/2) %*% bx[av, ]), S[av, av]^(1/2) %*% by[av])
    c(sqrt(t(mod$residuals) %*% (mod$residuals) / (mod$df.residual)), length(av))
  })
  rse_inc = rse[1, 2:lamlen] - rse[1, 1:(lamlen-1)]
  het = which(rse[1, 2:lamlen] > 1 & rse_inc > ((qchisq(0.95, 1) / rse[2, 2:lamlen])))
  if (length(het) == 0){
    lam_pos = 1
  } else {
    lam_pos = min(het)
  }
  num_valid = rev(sapply(1:lamlen, function(j){sum(alas$beta[, j]==0)}))
  min_lam_pos = min(which(num_valid > k))
  if (lam_pos < min_lam_pos){lam_pos = min_lam_pos}
  return(list(fit = alas$beta[, (lamlen - lam_pos + 1)], lambda = lamseq[lam_pos]))
}

mvmr_lasso = function(bx, by, seby){
  p = dim(as.matrix(bx))[1]
  k = dim(as.matrix(bx))[2]
  S = diag(seby^-2)
  sn = sign(bx[, 1])
  bx_or = bx * sn
  by_or = by * sn
  cv.alas = cv.mvmr_lasso(bx_or, by_or, seby)
  a1 = cv.alas$fit
  e = by_or - a1
  thest = solve(t(bx_or) %*% S %*% bx_or, t(bx_or) %*% S %*% e)
  v = which(a1==0)
  mvmr_mod = mr_mvivw(mr_mvinput(bx = bx_or[v, ], bxse = bx_or[v, ],
                                 by = by_or[v], byse = seby[v]))
  th_post = mvmr_mod$Estimate
  se_post = mvmr_mod$StdError
  pval<-p_fun(th_post,se_post)
  return(list(thest = thest, a = a1, lambda = cv.alas$lambda,
              th_post = th_post, se_post = se_post,pval=pval))
}
#mvmr_lasso(bx=bx_matrix,  by=beta_se_all$Stroke_b, seby=beta_se_all$Stroke_se)

mvmr_robust = function(bx, by, seby, k.max = 500, maxit.scale = 500){
  robmod = lmrob(by ~ bx - 1, weights = seby^-2, k.max = k.max,
                 maxit.scale = maxit.scale)
  coefficients = summary(robmod)$coef[, 1]
  se = summary(robmod)$coef[, 2] / min(summary(robmod)$sigma, 1)
  pval<-p_fun(coefficients,se)
  return(list("coefficients" = coefficients, "se" = se,"pval"=pval))
}
#mvmr_robust(bx=bx_matrix,  by=beta_se_all$Stroke_b, seby=beta_se_all$Stroke_se)

