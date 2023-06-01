simdata <- function(n,n_accrual, n_follow, rr, lambda2_diff, lambda1_diff) {

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

### output PBIR ###

permutation <- order(data$PFS)
perm <- data[permutation,]
perm <- perm[perm$censoring==1,]
p95 <- as.numeric(quantile(perm$PFS, 0.95))
fit1 <- mduration(data$PFS[data$trt==1], data$censoring[data$trt==1], data$TTR[data$trt==1], data$resp[data$trt==1], time.max = p95)
tmax1 <- fit1$time.truncation
fit0 <- mduration(data$PFS[data$trt==0], data$censoring[data$trt==0], data$TTR[data$trt==0], data$resp[data$trt==0], time.max=p95)
tmax0 <- fit0$time.truncation
diff=(fit1$meandor.est-fit0$meandor.est)
se=sqrt(fit1$meandor.se^2+fit0$meandor.se^2)
var = se^2
z=diff/se
ci=diff+c(-1,1)*qnorm(0.975)*se
pval_pbir=1-pchisq(z^2, 1)

### output tir ###

data$tir <- ifelse(data$resp==1, data$DOR, 0)
data$newcensor <- ifelse(data$resp==1,data$censoring,1)
tir_logrank <- survdiff(Surv(tir,newcensor)~trt, data=data)
cox_tir <-coxph(Surv(tir,newcensor)~trt,data=data)
pval_tir <- summary(cox_tir)$sctest[3]
coef_tir <- cox_tir$coef

### output dor ###

dor_logrank <- survdiff(Surv(DOR,censoring)~trt, data=data[data$resp==1,])
cox_dor <- coxph(Surv(DOR,censoring)~trt, data=data[data$resp==1,])
pval_dor <- summary(cox_dor)$sctest[3]
coef_dor <- cox_dor$coef

### output fisher ###

ctrl_repyes <- nrow(dt_c[dt_c$resp==1,])
ctrl_repno <- nrow(dt_c[dt_c$resp==0,])
exp_repyes <- nrow(dt_e[dt_e$resp==1,])
exp_repno <- nrow(dt_e[dt_e$resp==0,])
fishdat <- data.frame("repno"=c(ctrl_repno, exp_repno),
                      "repyes"=c(ctrl_repyes,exp_repyes),
                      row.names=c("Ctrl","Expé"),
                      stringsAsFactors = F)
test <- fisher.test(fishdat)
fish <- log(test$estimate)
pval_fish <- test$p.value

### output final ###

result <- cbind(cens_rate, diff, var, pval_pbir, coef_tir, pval_tir, coef_dor, pval_dor, fish, pval_fish, tmax1,tmax0)
print(round(result, 4))
}

# output final

diff <- matrix(NA,11,11)
cens <- matrix(NA,11,11)
vari <- matrix(NA,11,11)
pval_pbir <- matrix(NA,11,11)
coef_tir <- matrix(NA,11,11)
pval_tir <- matrix(NA,11,11)
coef_dor <- matrix(NA,11,11)
pval_dor <- matrix(NA,11,11)
fish <- matrix(NA,11,11)
pval_fish <- matrix(NA,11,11)
tmax1 <- matrix(NA,11,11)
tmax0 <- matrix(NA,11,11)

rr_var <- c(-0.25,-0.20, -0.15, -0.1,-0.05,0,0.05,0.1,0.15,0.2,0.25)
lambda2_var <- c(1.5,1.4,1.3,1.2,1.1,1,0.9,0.8,0.7,0.6,0.5)
for (i in 1:11) {
  for (j in 1:11) {
    set.seed(123)
    B <- 1000
    simlist <- (data.frame(t(replicate(n=B, simdata(260,20,7,rr_var[i],lambda2_var[j],1),simplify=T))))
    cens[j,i] <- mean(simlist[,1])
    diff[j,i] <- mean(simlist[,2])
    vari[j,i] <- mean(simlist[,3])
    pval_pbir[j,i] <- sum(simlist[,4] < 0.05)/B
    coef_tir[j,i] <- exp(mean(simlist[,5]))
    pval_tir[j,i] <- sum(simlist[,6] < 0.05)/B
    coef_dor[j,i] <- exp(mean(simlist[,7]))
    pval_dor[j,i] <- sum(simlist[,8] < 0.05)/B
    fish[j,i] <- exp(mean(simlist[,9]))
    pval_fish[j,i] <- sum(simlist[,10] < 0.5)/B
    tmax1[j,i] <- mean(simlist[,11])
    tmax0[j,i] <- mean(simlist[,12])
  }
}

colnames(cens) <- rr_var
rownames(cens) <- lambda2_var
cens

colnames(diff) <- rr_var
rownames(diff) <- lambda2_var
diff

colnames(vari) <- rr_var
rownames(vari) <- lambda2_var
vari

colnames(pval_pbir) <- rr_var
rownames(pval_pbir) <- lambda2_var
pval_pbir

colnames(coef_tir) <- rr_var                                                                                                                                                                                                                                                                                                                                                   
rownames(coef_tir) <- lambda2_var
coef_tir

colnames(pval_tir) <- rr_var
rownames(pval_tir) <- lambda2_var
pval_tir

colnames(coef_dor) <- rr_var
rownames(coef_dor) <- lambda2_var
coef_dor

colnames(pval_dor) <- rr_var
rownames(pval_dor) <- lambda2_var
pval_dor

colnames(fish) <- rr_var
rownames(fish) <- lambda2_var
fish

colnames(pval_fish) <- rr_var
rownames(pval_fish) <- lambda2_var
pval_fish


# ggplot

theme_set(theme_bw())

diff2 <- data.frame(diff)
colnames(diff2) <- rr_var
diff2 <- cbind(lambda2=lambda2_var,diff2)
diff_tt <- pivot_longer(diff2, !lambda2, names_to='p',values_to='diff_in_PBIR')
ggplot(data=diff_tt,aes(x=lambda2,y=diff_in_PBIR, color=p)) + geom_line() + scale_color_discrete(breaks=c(-0.25, -0.20, -0.15, -0.1,-0.05,0,0.05,0.1,0.15,0.2,0.25)) + labs(color = TeX(r'($p_{diff}$)'), y= "Difference in PBIR (in days)", x = TeX(r'($\lambda_{2diff}$)'), title='Difference in PBIR between experimental and control groups',subtitle=TeX(r'(With different values of $p_{diff}$ and $\lambda_{2diff}$)'))

coef_dor2 <- data.frame(coef_dor)
colnames(coef_dor2) <- rr_var
coef_dor2 <- cbind(lambda2=lambda2_var,coef_dor2)
coef_dor_tt <- pivot_longer(coef_dor2, !lambda2, names_to='p',values_to='HR_DOR')
ggplot(data=coef_dor_tt,aes(x=lambda2,y=HR_DOR, color=p)) + geom_line() + scale_color_discrete(breaks=c(-0.25, -0.20, -0.15, -0.1,-0.05,0,0.05,0.1,0.15,0.2,0.25)) + labs(color = TeX(r'($p_{diff}$)'), y= "HR of DOR", x = TeX(r'($lambda_{2diff}$)'), title='Hazard Ratio of DOR between experimental and control groups',subtitle=TeX(r'(With different values of $p_{diff}$ and $\lambda_{2diff}$)'))

coef_tir2 <- data.frame(coef_tir)
colnames(coef_tir2) <- rr_var
coef_tir2 <- cbind(lambda2=lambda2_var,coef_tir2)
coef_tir_tt <- pivot_longer(coef_tir2, !lambda2, names_to='p',values_to='HR_TIR')
ggplot(data=coef_tir_tt,aes(x=lambda2,y=HR_TIR, color=p)) + geom_line() + scale_color_discrete(breaks=c(-0.25, -0.20, -0.15, -0.1,-0.05,0,0.05,0.1,0.15,0.2,0.25)) + labs(color = TeX(r'($p_{diff}$)'), y= "HR of TIR", x = TeX(r'($\lambda_{2diff}$)'),title='Hazard Ratio of TIR between experimental and control groups',subtitle=TeX(r'(With different values of $p_{diff}$ and $\lambda_{2diff}$)'))

### contour ###

diff_a <- as.matrix(arrange(diff2,lambda2))
rownames(diff_a) <- diff_a[,1]
diff_a <- diff_a[,-1]
contour(x=as.numeric(rownames(diff_a)),y=as.numeric(colnames(diff_a)),diff_a, nlevels=10, plot.title = title(main = "Difference in PBIR (days) between experimental and control groups", xlab=TeX(r'($\lambda_{2diff}$)'),  ylab=TeX(r'($p_{diff}$)')))
abline(v=c(1,0.6,0.8,1.2,1.4),col=c("blue","grey", "grey","grey", "grey"),lty=2)
abline(h=c(0,-0.2,-0.1,0.1,0.2),col=c("blue","grey", "grey","grey", "grey"),lty=2)
filled.contour(x=as.numeric(colnames(diff_a)),y=as.numeric(rownames(diff_a)),diff_a, nlevels=10, plot.title = title(main = "Difference in PBIR (days) between experimental and control groups", xlab=TeX(r'($p_{diff}$)'),  ylab=TeX(r'($\lambda_{2diff}$)')))

### suite ggplot

cens2 <- data.frame(cens)
colnames(cens2) <- rr_var
cens2 <- cbind(lambda2=lambda2_var,cens2)
cens_tt <- pivot_longer(cens2, !lambda2, names_to='p',values_to='censoring_rate')
ggplot(data=cens_tt,aes(x=lambda2,y=censoring_rate, color=p)) + geom_line() + scale_color_discrete(breaks=c(-0.25, -0.20, -0.15, -0.1,-0.05,0,0.05,0.1,0.15,0.2,0.25)) + labs(color = TeX(r'($p_{diff}$)'), y= "Censoring rate", x = TeX(r'($\lambda_{2diff}$)'), title='Censoring rate',subtitle=TeX(r'(With different values of $p_{diff}$ and $\lambda_{2diff}$)'))

vari2 <- data.frame(vari)
colnames(vari2) <- rr_var
vari2 <- cbind(lambda2=lambda2_var,vari2)
vari_tt <- pivot_longer(vari2, !lambda2, names_to='p',values_to='variance')
ggplot(data=vari_tt,aes(x=lambda2,y=variance, color=p)) + geom_line() + scale_color_discrete(breaks=c(-0.25, -0.20, -0.15, -0.1,-0.05,0,0.05,0.1,0.15,0.2,0.25)) + labs(color = TeX(r'($p_{diff}$)'), y= "Variance of difference in PBIR (days)", x = TeX(r'($\lambda_{2diff}$)'), title='Variance of difference in PBIR between experimental and control groups',subtitle=TeX(r'(With different values of $p_{diff}$ and $\lambda_{2diff}$)'))

fish2 <- data.frame(fish)
colnames(fish2) <- rr_var
fish2 <- cbind(lambda2=lambda2_var,fish2)
fish_tt <- pivot_longer(fish2, !lambda2, names_to='p',values_to='OR_fisher')
ggplot(data=fish_tt,aes(x=lambda2,y=OR_fisher, color=p)) + geom_line() + scale_color_discrete(breaks=c(-0.25, -0.20, -0.15, -0.1,-0.05,0,0.05,0.1,0.15,0.2,0.25)) + labs(color = TeX(r'($p_{diff}$)'), y= "OR", x = TeX(r'($\lambda_{2diff}$)'), title='OR between experimental and control groups',subtitle=TeX(r'(With different values of $p_{diff}$ and $\lambda_{2diff}$)'))

pval_pbirb <- data.frame(pval_pbir)
colnames(pval_pbirb) <- rr_var
pval_pbirb <- cbind(lambda2=lambda2_var,pval_pbirb)
pval_pbir_tt <- pivot_longer(pval_pbirb, !lambda2, names_to='p',values_to='pval_pbir')
ggplot(data=pval_pbir_tt,aes(y=pval_pbir,x=lambda2, color=p)) + geom_line() + scale_color_discrete(breaks=c(-0.25, -0.20, -0.15, -0.1,-0.05,0,0.05,0.1,0.15,0.2,0.25)) + labs(color = TeX(r'($p_{diff}$)'), y= "Pr(pval_pbir < 0.05)", x = TeX(r'($\lambda_{2diff}$)'), title='Proportion of trials with significant p-value for PBIR method',subtitle=TeX(r'(With different values of $p_{diff}$ and $\lambda_{2diff}$)'))

pval_tir2 <- data.frame(pval_tir)
colnames(pval_tir2) <- rr_var
pval_tir2 <- cbind(lambda2=lambda2_var,pval_tir2)
pval_tir_tt <- pivot_longer(pval_tir2, !lambda2, names_to='p',values_to='pval_tir')
ggplot(data=pval_tir_tt,aes(y=pval_tir,x=lambda2, color=p)) + geom_line() + scale_color_discrete(breaks=c(-0.25, -0.20, -0.15, -0.1,-0.05,0,0.05,0.1,0.15,0.2,0.25)) + labs(color = TeX(r'($p_{diff}$)'), y= "Pr(pval_tir < 0.05)", x = TeX(r'($\lambda_{2diff}$)'), title='Proportion of trials with significant p-value for TIR method',subtitle=TeX(r'(With different values of $p_{diff}$ and $\lambda_{2diff}$)'))

pval_dor2 <- data.frame(pval_dor)
colnames(pval_dor2) <- rr_var
pval_dor2 <- cbind(lambda2=lambda2_var,pval_dor2)
pval_dor_tt <- pivot_longer(pval_dor2, !lambda2, names_to='p',values_to='pval_dor')
ggplot(data=pval_dor_tt,aes(y=pval_dor,x=lambda2, color=p)) + geom_line() + scale_color_discrete(breaks=c(-0.25, -0.20, -0.15, -0.1,-0.05,0,0.05,0.1,0.15,0.2,0.25)) + labs(color = TeX(r'($p_{diff}$)'), y= "Pr(pval_dor < 0.05)", x = TeX(r'($\lambda_{2diff}$)'), title='Proportion of trials with significant p-value for DOR method',subtitle=TeX(r'(With different values of $p_{diff}$ and $\lambda_{2diff}$)'))

pval_fish2 <- data.frame(pval_fish)
colnames(pval_fish2) <- rr_var
pval_fish2 <- cbind(lambda2=lambda2_var,pval_fish2)
pval_fish_tt <- pivot_longer(pval_fish2, !lambda2, names_to='p',values_to='pval_fish')
ggplot(data=pval_fish_tt,aes(y=pval_fish,x=lambda2, color=p)) + geom_line() + scale_color_discrete(breaks=c(-0.25, -0.20, -0.15, -0.1,-0.05,0,0.05,0.1,0.15,0.2,0.25))+ labs(color = TeX(r'($p_{diff}$)'), y= "Pr(pval_fisher < 0.05)", x = TeX(r'($\lambda_{2diff}$)'), title='Proportion of trials with significant p-value for Fisher Test',subtitle=TeX(r'(With different values of $p_{diff}$ and $\lambda_{2diff}$)'))

# export table

library(tidyverse)
write.csv(rev(data.frame(t(round(diff,2)))), "diff.csv")
write.csv(rev(data.frame(t(round(cens,4)))), "cens.csv")
write.csv(rev(data.frame(t(vari))), "vari.csv")
write.csv(rev(data.frame(t(round(pval_pbir,4)))), "pval_pbir.csv")
write.csv(rev(data.frame(t(round(coef_tir,2)))), "coef_tir.csv")
write.csv(rev(data.frame(t(round(pval_tir,4)))), "pval_tir.csv")
write.csv(rev(data.frame(t(round(coef_dor,2)))), "coef_dor.csv")
write.csv(rev(data.frame(t(round(pval_dor,4)))), "pval_dor.csv")
write.csv(rev(data.frame(t(round(fish,4)))), "fish.csv")
write.csv(rev(data.frame(t(round(pval_fish,4)))), "pval_fish.csv")

#graph pval methods

p_p <- pval_pbir[,9]
p_d <- pval_dor[,9]
p_t <- pval_tir[,9]
p_curve <- cbind(p_p,p_d,p_t)
p_curve2 <- data.frame(p_curve)
p_curve2b<- cbind(lambda2=lambda2_var,p_curve2)
p_curve_tt <- pivot_longer(p_curve2b, !lambda2,names_to='Method',values_to='pval')
p_curve_tt$Method <- factor(p_curve_tt$Method, levels = c("p_d", "p_t", "p_p"))
ggplot(data=p_curve_tt,aes(y=pval,x=lambda2, color=Method)) + geom_line()+ labs(color = "Method", y= "Pr(pval_pbir < 0.05)", x = TeX(r'($\lambda_{2diff}$)'), title='Proportion of trials with significant p-value for three methods',subtitle=TeX(r'(With value of $p_{diff}=0.15$)')) + scale_color_discrete(name = "Method", labels = c("DOR", "TIR","PBIR"))

