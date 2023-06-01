simdata_lambda1 <- function(n,n_accrual, n_follow, rr, lambda2_diff, lambda1_diff) {
  
  ##### CTRL ARM ####
  
  ctrl <- aeendpt[aeendpt$TRTCD==2,]
  id <- 1:n
  responders <- ctrl[ctrl$ROSLR==1,]
  nresponders <- ctrl[ctrl$ROSLR==0,]
  p_c <- nrow(responders)/260
  resp <- rbinom(n,1,p_c)
  dt_c <- data.frame(id,resp)
  dt_c$trt <- 0
  respn <- nrow(dt_c[dt_c$resp==1,])
  TTR_NEM <- responders$RDYLR
  TTR <- sample(x=TTR_NEM,size=respn, replace=T)
  dt_c$TTR[dt_c$resp==1] <- TTR
  dt_c$TTR[dt_c$resp==0] <- NA
  
  #lambda 2#
  
  lambda2_fit <- survfit(Surv(RDURLR,RDURLRC)~1,data=ctrl)
  lambda2_c <- sum(lambda2_fit$n.event)/sum(lambda2_fit$time)
  DOR <- rexp(respn,rate=lambda2_c)
  dt_c$DOR[dt_c$resp==1] <- DOR
  dt_c$DOR[dt_c$resp==0] <- NA
  
  #lamba 1
  
  lambda1_fit <- survfit(Surv(PFSDYLR,PFSLR)~1,data=nresponders)
  lambda1_c <- sum(lambda1_fit$n.event)/sum(lambda1_fit$time)
  PFS <- rexp(n-respn,rate=lambda1_c)
  dt_c$prog[dt_c$resp==1] <- DOR+TTR
  dt_c$prog[dt_c$resp==0] <- PFS
  
  #accrual
  
  dt_c$accrual <- runif(n,0,n_accrual*30.4375)
  dt_c$Time <- dt_c$prog + dt_c$accrual
  
  #censoring
  
  dt_c$censoring <- ifelse(dt_c$Time<(n_accrual+n_follow)*30.4375, 1, 0)
  dt_c$PFS <- ifelse(dt_c$censoring==1, dt_c$prog, ((n_accrual+n_follow)*30.4375)-dt_c$accrual)
  
  
  #### EXPERIMENTAL ARM ####
  
  p_e <- p_c+rr 
  resp <- rbinom(n,1,p_e)
  dt_e <- data.frame(id, resp)
  dt_e$trt <- 1
  
  respn_e <- nrow(dt_e[dt_e$resp==1,])
  TTR_e <- sample(x=TTR_NEM,size=respn_e, replace=T)
  dt_e$TTR[dt_e$resp==1] <- TTR_e
  dt_e$TTR[dt_e$resp==0] <- NA
  
  #lambda 2
  
  lambda2_e <- lambda2_diff*lambda2_c
  DOR_e <- rexp(respn_e,rate=lambda2_e)
  dt_e$DOR[dt_e$resp==1] <- DOR_e
  dt_e$DOR[dt_e$resp==0] <- NA
  
  #lamdba 1
  
  lambda1_e <- lambda1_diff*lambda1_c
  PFS_e <- rexp(n-respn_e,rate=lambda1_e)
  dt_e$prog[dt_e$resp==1] <- DOR_e+TTR_e
  dt_e$prog[dt_e$resp==0] <- PFS_e
  
  #accrual
  
  dt_e$accrual <- runif(n,0,n_accrual*30.4375)
  dt_e$Time <- dt_e$prog + dt_e$accrual
  
  #censoring
  
  dt_e$censoring <- ifelse(dt_e$Time<(n_accrual+n_follow)*30.4375, 1, 0)
  dt_e$PFS <- ifelse(dt_e$censoring==1, dt_e$prog, ((n_accrual+n_follow)*30.4375)-dt_e$accrual)
  
  #### MERGE DATASETS ####
  
  data <- rbind(dt_c,dt_e)
  
  ### output censoring ###
  
  cens_rate <- nrow(data[data$censoring==0,])/nrow(data)
  exp <- data[data$trt==1,]
  exp_nresp <- exp[exp$resp==0,]
  cens_expnresp <- nrow(exp_nresp[exp_nresp$censoring==0,])/nrow(exp_nresp)
  nexp <- data[data$trt==0,]
  nexp_nresp <- nexp[nexp$resp==0,]
  cens_nexpnresp <- nrow(nexp_nresp[nexp_nresp$censoring==0,])/nrow(nexp_nresp)
  
  ### output PBIR 95_resp ###
  
  permutation <- order(data$PFS)
  perm <- data[permutation,]
  perm <- perm[perm$censoring==1,]
  perm_resp <- perm[perm$resp==1,]
  p95_resp <- as.numeric(quantile(perm_resp$PFS, 0.95))
  fit95_resp_1 <- mduration(data$PFS[data$trt==1], data$censoring[data$trt==1], data$TTR[data$trt==1], data$resp[data$trt==1],time.max = p95_resp)
  fit95_resp_0 <- mduration(data$PFS[data$trt==0], data$censoring[data$trt==0], data$TTR[data$trt==0], data$resp[data$trt==0],time.max= p95_resp)
  diff95_resp=(fit95_resp_1$meandor.est-fit95_resp_0$meandor.est)
  tmax95_resp <- fit95_resp_1$time.truncation
  
  ### output dor ###
  
  dor_logrank <- survdiff(Surv(DOR,censoring)~trt, data=data[data$resp==1,])
  cox_dor <- coxph(Surv(DOR,censoring)~trt, data=data[data$resp==1,])
  dor <- cox_dor$coef
  
  ### output tir ###
  
  data$tir <- ifelse(data$resp==1, data$DOR, 0)
  data$newcensor <- ifelse(data$resp==1,data$censoring,1)
  tir_logrank <- survdiff(Surv(tir,newcensor)~trt, data=data)
  cox_tir <-coxph(Surv(tir,newcensor)~trt,data=data)
  tir <- cox_tir$coef
  
  ### output final ###
  
  result <- cbind(cens_rate, diff95_resp, tmax95_resp, dor, tir)
  print(round(result, 4))
  
  #data_lambda1 <<- data
}


lambda1_var <- c(1.5,1.4,1.3,1.2,1.1,1,0.9,0.8,0.7,0.6,0.5)
mfu_var <- c(1,3,5,7,9,11,13,100)
cens_lambda1 <- matrix(NA,11,8)
diff95_resp_lambda1 <- matrix(NA,11,8)
tmax95_resp_lambda1 <- matrix(NA,11,8)
coef_dor_lambda1 <- matrix(NA,11,8)
coef_tir_lambda1 <- matrix(NA,11,8)

for (i in 1:8) {
  for (j in 1:11) {
    set.seed(123)
    B <- 1000
    simlist <- (data.frame(t(replicate(n=B, simdata_lambda1(260,20,mfu_var[i],0.15,0.8,lambda1_var[j]),simplify=T))))
    cens_lambda1[j,i] <- mean(simlist[,1])
    diff95_resp_lambda1[j,i] <- mean(simlist[,2])
    tmax95_resp_lambda1[j,i] <- mean(simlist[,3])
    coef_dor_lambda1[j,i] <- exp(mean(simlist[,4]))
    coef_tir_lambda1[j,i] <- exp(mean(simlist[,5]))
  }
}

diff95_resp_2_lambda1 <- data.frame(diff95_resp_lambda1)
colnames(diff95_resp_2_lambda1) <- mfu_var
diff95_resp_2_lambda1 <- cbind(lambda1=lambda1_var,diff95_resp_2_lambda1)
diff95_resp_tt_lambda1 <- pivot_longer(diff95_resp_2_lambda1, !lambda1, names_to='mfu',values_to='diff_in_PBIR')
ggplot(data=diff95_resp_tt_lambda1,aes(x=lambda1,y=diff_in_PBIR, color=mfu)) + geom_line() + scale_color_discrete(breaks=c(1,3,5,7,9,11,13,100)) + labs(color = TeX(r'($mFU (months)$)'), y= "Difference in PBIR (in days)", x = TeX(r'($\lambda_{1diff}$)'), title='Difference in PBIR between experimental and control groups',subtitle=TeX(r'(With different values of $mFU$ and $\lambda_{1diff}$)'))

dor2_lambda1 <- data.frame(coef_dor_lambda1[,4])
dor2_lambda1 <- cbind(lambda1=lambda1_var,dor2_lambda1)
dor_tt_lambda1 <- pivot_longer(dor2_lambda1, !lambda1, values_to='HR_DOR')
ggplot(data=dor_tt_lambda1,aes(y=HR_DOR,x=lambda1)) + geom_line(color='#C77CFF') + labs(y= "HR of DOR", x = TeX(r'($\lambda_{1diff}$)'), title='Hazard Ratio of DOR between experimental and control groups',subtitle=TeX(r'(With different values of $\lambda_{1diff}$)'))

tir2_lambda1 <- data.frame(coef_tir_lambda1[,4])
tir2_lambda1 <- cbind(lambda1=lambda1_var,tir2_lambda1)
tir_tt_lambda1 <- pivot_longer(tir2_lambda1, !lambda1, values_to='HR_TIR')
ggplot(data=tir_tt_lambda1,aes(y=HR_TIR,x=lambda1)) + geom_line(color='#C77CFF') + labs(y= "HR of TIR", x = TeX(r'($\lambda_{1diff}$)'), title='Hazard Ratio of TIR between experimental and control groups',subtitle=TeX(r'(With different values of $\lambda_{1diff}$)'))

rev(data.frame(t(round(diff95_resp_lambda1,4))))
