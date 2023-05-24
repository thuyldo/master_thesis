library(survival)
library(PBIR)
library(purrr)
library(ggplot2)
library(tidyr)
library(latex2exp)
library(dplyr)

dataset <- function(n,n_accrual, n_follow, rr, lambda2_diff, lambda1_diff) {
  
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
  
  #lambda 2
  
  lambda2_fit <- survfit(Surv(RDURLR,RDURLRC)~1,data=ctrl)
  lambda2_c <- sum(lambda2_fit$n.event)/sum(lambda2_fit$time)
  DOR <- rexp(respn,rate=lambda2_c)
  dt_c$DOR[dt_c$resp==1] <- DOR
  dt_c$DOR[dt_c$resp==0] <- NA
  
  #lambda 1
  
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
  
  #lambda 2#
  
  lambda2_e <- lambda2_diff*lambda2_c
  DOR_e <- rexp(respn_e,rate=lambda2_e)
  dt_e$DOR[dt_e$resp==1] <- DOR_e
  dt_e$DOR[dt_e$resp==0] <- NA
  
  #lamba 1
  
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
  
  # MERGE DATASETS ##
  
  data <- rbind(dt_c,dt_e)
  
  data_pbir <<- data
}

set.seed(123)

data.frame(t(replicate(n=10, dataset(260,20,7,0.2,0.7,1),simplify=T)))

# Experimental arm : PD/PDR curves

pbir_expe <- data_pbir[data_pbir$trt==1,]
PFS <- survfit(Surv(pbir_expe$PFS/30.4375,pbir_expe$censoring)~1)

pbir_expe$PD <- pmin(pbir_expe$PFS,pbir_expe$TTR,na.rm=T)
pbir_expe$newcensor <- 1*(pbir_expe$resp+pbir_expe$censoring>0)
PD <- survfit(Surv(pbir_expe$PD/30.4375,pbir_expe$newcensor)~1)
plot(PD,xlab="time (months)",  ylab="Survival",main="P/D-free and P/D/R-free curves", col='#009999', xlim=c(0,20), conf.int = F)
lines(PFS, col="#9999FF", conf.int = FALSE)
legend(x=10,y=1, legend=c("P/D-free", "P/D/R-free"),fill=c("#9999FF", "#009999"))

polygon(c(PFS$time,max(PFS$time),min(PFS$time)),c(PFS$surv,0,0),col="#CCCCCC", density=NA)
polygon(c(PD$time,max(PD$time),min(PD$time)),c(PD$surv,0,0),col="white",density=NA)

lines(PFS, col="#9999FF", conf.int = FALSE)
lines(PD, col="#009999", conf.int = FALSE)

# Expe arm : PBIR

PBIRfit <- PBIR1(pbir_expe$PFS/30.4375,pbir_expe$censoring,pbir_expe$TTR/30.4375,pbir_expe$resp, time=NULL,alpha=0.95)
tt1=c(0, PBIRfit$time)
PBIR=c(0, PBIRfit$PBIR)
B1=length(tt1)
tt1=rep(tt1, rep(2, B1))[-1]
PBIR=rep(PBIR, rep(2, B1))[-(2*B1)]
plot(range(PBIRfit$time), range(PBIRfit$PBIR),
     xlab="time (months)",  ylab="PBIR",
     main="Probability of Being in Response", type="n")
lines(tt1, PBIR, col="blue")


PPPP <- cbind(PFS$time,PFS$surv)
colnames(PPPP) <- c("time","PFS")
PDDD <- cbind(PD$time, PD$surv)
colnames(PDDD) <- c("time", "PDR")
PPDD <- merge(PPPP,PDDD, by="time")
PPDD$PFS_PDR <- PPDD$PFS - PPDD$PDR
PBBB <- cbind(PBIRfit$time, PBIRfit$PBIR)
colnames(PBBB) <- c("time","PBIR")
PPDD <- merge(PPDD,PBBB,by="time")

# 2 groups comparison
PBIR_pbir <- function(T2,deltaT2,T1,deltaT1,group1, group2){
  PBIRfit1 <- PBIR1(T2[group1],deltaT2[group1],T1[group1],deltaT1[group1], time=NULL,alpha=0.95)
  PBIRfit0 <- PBIR1(T2[group2],deltaT2[group2],T1[group2],deltaT1[group2], time=NULL,alpha=0.95)
  tt1=c(0, PBIRfit1$time)
  PBIR1=c(0, PBIRfit1$PBIR)
  B1=length(tt1)
  tt1=rep(tt1, rep(2, B1))[-1]
  PBIR1=rep(PBIR1, rep(2, B1))[-(2*B1)]
  tt0=c(0, PBIRfit0$time)
  PBIR0=c(0, PBIRfit0$PBIR)
  B0=length(tt0)
  tt0=rep(tt0, rep(2, B0))[-1]
  PBIR0=rep(PBIR0, rep(2, B0))[-(2*B0)]
  plot(range(c(PBIRfit1$time, PBIRfit0$time)), range(c(PBIRfit1$PBIR, PBIRfit0$PBIR)),
       xlab="time (months)",  ylab="PBIR",
       main="Probability of Being in Response", type="n")
  lines(tt0, PBIR0, col="red",lty="dashed")
  lines(tt1, PBIR1, col="blue")
  legend(x=11,y=0.4, legend=c("Control arm","Experimental arm"),fill=c("red","blue"), bty="n")}
PBIR_pbir(data_pbir$PFS/30.4375,data_pbir$censoring,data_pbir$TTR/30.4375,data_pbir$resp,data_pbir$trt==1,data_pbir$trt==0)

# diff

fit2=PBIR2(data_pbir$PFS/30.4375, data_pbir$censoring, data_pbir$TTR/30.4375, data_pbir$resp, data_pbir$trt)

tt=fit2$time
diff=fit2$diff
low=fit2$ci.low
up=fit2$ci.up
B=length(tt)+1
tt=c(0, tt)
diff=c(0, diff)
low=c(0, low)
up=c(0, up)
tt=rep(tt, rep(2, B))[-1]
diff=rep(diff, rep(2, B))[-(2*B)]
low=rep(low, rep(2, B))[-(2*B)]
up=rep(up, rep(2, B))[-(2*B)]
plot(range(c(tt, 0)), range(c(low, up)), xlab="time (months)", ylab="Difference in PBIR", main='Difference in PBIR with 95% Confidence Interval', lwd=2, type="n", xlim=c(0,20))
lines(tt, diff, lwd=2, col=3)
lines(tt, low, col=2)
lines(tt, up, col=2)
lines(range(fit2$time), rep(0, 2), col=4, lty=4)
