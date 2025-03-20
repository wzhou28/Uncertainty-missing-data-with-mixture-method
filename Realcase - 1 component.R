library(ssmrob)
library(ars) #for ARS
library(pracma) #for the "erf" function
library(doParallel)
library (plyr)
library(data.table) 
library(sampleSelection)

data(MEPS2001)
v <- 500  #ARS sampling N
registerDoParallel(cores=4)
###################################################################################################################################
test_1 <- MEPS2001

test_1$R <- 1
test_1$gender <- 0
test_1$race <- 0
test_1$insu <- 0

for (i in 1:length(test_1$educ)){
  if (test_1$ambexp[i]==0) test_1$R[i] <- 0
  if (test_1$female[i]=='TRUE') test_1$gender[i] <- 1
  if (test_1$blhisp[i]=='TRUE') test_1$race[i] <- 1
  if (test_1$ins[i]=='TRUE') test_1$insu[i] <- 1
}
test_1 <- subset(test_1, income>0) #delete 5 observations with <=0 income (since we will log it later)
# test_1$income <- log(test_1$income)
colnames(test_1)[8] <- c("exp")
colnames(test_1)[5] <- c("disease")
test_1 <- test_1[,-c(4,6,7,9,10,11,12)]
#age: age is in 10-year
#female: 0/1
#educ: 
#blhisp: 0/1

#heckman
# summary(selection(R~XS1+XS2, Y~XO1+XO2, data=test))
# summary(selection(R~income+age+gender+educ+race+disease+insu, 
#                   exp~     age+gender+educ+race+disease+insu, data=test_1))
# test_1 <- MEPS2001
# test_1$R <- 1
# for (i in 1:length(test_1$educ)){
#   if (test_1$ambexp[i]==0) test_1$R[i] <- 0
# }
# summary(selection(R~income+age+female+educ+blhisp+totchr+ins, 
#                   lambexp~     age+female+educ+blhisp+totchr+ins, data=test_1))

######################################################################################################
# hist(test_1$age)
# table(test_1$gender)
# table(test_1$race)
# hist(test_1$educ)
# hist((test_1$income))
# hist(log(test_1$income))
# table(test_1$insu)
# table(test_1$disease)
# hist(test_1$exp)

para <- data.frame(iteration=1:2000,
                   # prop, #for binary cluster var "pik" (only consider two clusters)
                   #all below parameters are paired (for two clusters)
                   u_a=NA, sd_a=NA, #age
                   u_g=NA, #gender
                   u_r=NA, #race
                   u_e=NA, sd_e=NA,  #edu
                   b0_inc=NA, b1_inc=NA, b2_inc=NA, b3_inc=NA, b4_inc=NA, sigma_inc=NA,  #for income~age+gender+race+educ
                   b0_ins=NA, b1_ins=NA, #for insu~income
                   b0_d=NA, b1_d=NA, b2_d=NA, b3_d=NA, b4_d=NA, #disease~age+gender+race+insu
                   b0_e=NA, b1_e=NA, b2_e=NA, b3_e=NA, b4_e=NA, b5_e=NA, b6_e=NA, sigma_e=NA,  #for exp~age+gender+educ+race+disease+insu
                   a0_e=NA, a1_e=NA, a2_e=NA, a3_e=NA, a4_e=NA, a5_e=NA, a6_e=NA, a7_e=NA, a8_e=NA) #R~income+age+gender+educ+race+disease+insu

c1 <- coefficients(glm(formula=log(income)~age+gender+race+educ, data=test_1)) #note, income~logNormal
c2 <- coefficients(glm(formula=insu~income, data=test_1, family = binomial(link = "logit")))
c3 <- coefficients(glm(formula=disease~age+gender+race+insu, data=test_1,family = poisson))
c4 <- coefficients(glm(formula=exp~age+gender+educ+race+disease+insu, data=test_1))
c5 <- c(-0.100, 0.003, 0.113, 0.700, 0.066, -0.400, 0.860, 0.166, -0.120) #note, now we don't have full data to compute coeff, so we just use Miao's est
  
#for non-mixture, no iteration on below para except for exp~ and missing model (since no imputation for those var)
#i.e., only c4 and c5 need to be updated
#for mixture, all need to be updated
para[,8:12] <- c1
para[,13] <- sqrt(var(lm(formula=log(income)~age+gender+race+educ, test_1)$residuals))
para[,14:15] <- c2
para[,16:20] <- c3
para[1,21:27] <- c4 
para[1,28] <- sqrt(var(lm(formula=exp~age+gender+educ+race+disease+insu, test_1)$residuals))
para[1,29:37] <- c5

para[,2] <- mean(test_1$age)
para[,3] <- sd(test_1$age)
para[,4] <- mean(test_1$gender)
para[,5] <- mean(test_1$race)
para[,6] <- mean(test_1$edu)
para[,7] <- sd(test_1$edu)
#note, although parameters for marginal Xobs part is directly identifiable, it is NOT directly identifiable under mixture

n <- length(test_1$age)

#stopped at k=1512
#EM starts here:
#E-step
test_t <- test_1
for (k in 1511:1511) { #iteration for EM
  
  if (k %% 100 ==0) {print(paste("EM-iter", k, "is running"))}
  test_1 <- test_t[,1:9]
  
  #如果中断了k，重新跑的话要把最开始的test_1的部分也重新跑一下
  #E-step 
  res <- foreach(i=1:n) %dopar% {
    # for (i in 1:n) {
    
    #for those with missing exp, we need to sample exp
    if (test_1$R[i]==0){ #R=0 means missing
      
      #2) sample exp, we use ADS (since P(exp|age+gender+educ+race+disease+insu,R) prop to 
      #                                 P(exp|age+gender+educ+race+disease+insu  )*P(R|exp,income+age+gender+educ+race+disease+insu))
      #current EM iteration parameter value 
      b0_e <- para[k,21]
      b1_e <- para[k,22]
      b2_e <- para[k,23]
      b3_e <- para[k,24]
      b4_e <- para[k,25]
      b5_e <- para[k,26]
      b6_e <- para[k,27]
      sigma_e <- para[k,28]
      a0_e  <- para[k,29]
      a1_e  <- para[k,30] 
      a2_e  <- para[k,31] 
      a3_e  <- para[k,32] 
      a4_e  <- para[k,33] 
      a5_e  <- para[k,34] 
      a6_e  <- para[k,35] 
      a7_e  <- para[k,36]
      a8_e  <- para[k,37]
      
      #draw "500" sample for missing exp
      income <- test_1$income[i]
      age <- test_1$age[i]
      gender <- test_1$gender[i]
      educ <- test_1$educ[i]
      race <- test_1$race[i]
      disease <- test_1$disease[i]
      insu <- test_1$insu[i]
      r <- test_1$R[i]
      f <- function(x3, r, income, age, gender, educ, race, disease, insu, b0_e, b1_e, b2_e, b3_e, b4_e, b5_e, b6_e, sigma_e, a0_e, a1_e, a2_e, a3_e, a4_e, a5_e, a6_e, a7_e, a8_e){
        -(x3-b0_e-b1_e*age-b2_e*gender-b3_e*educ-b4_e*race-b5_e*disease-b6_e*insu)^2/(2*sigma_e^2) + 
          r*log(pnorm(a0_e+a1_e*income+a2_e*age+a3_e*gender+a4_e*educ+a5_e*race+a6_e*disease+a7_e*insu+a8_e*x3))+(1-r)*log(1-pnorm(a0_e+a1_e*income+a2_e*age+a3_e*gender+a4_e*educ+a5_e*race+a6_e*disease+a7_e*insu+a8_e*x3))} #log joint f of exp ~Normal
      fprima <- function(x3, r, income, age, gender, educ, race, disease, insu, b0_e, b1_e, b2_e, b3_e, b4_e, b5_e, b6_e, sigma_e, a0_e, a1_e, a2_e, a3_e, a4_e, a5_e, a6_e, a7_e, a8_e){
        -(x3-b0_e-b1_e*age-b2_e*gender-b3_e*educ-b4_e*race-b5_e*disease-b6_e*insu)/(sigma_e^2) + 
          sqrt(2/pi)*a8_e*r*exp(-(a0_e+a1_e*income+a2_e*age+a3_e*gender+a4_e*educ+a5_e*race+a6_e*disease+a7_e*insu+a8_e*x3)^2/2)/(erf((a0_e+a1_e*income+a2_e*age+a3_e*gender+a4_e*educ+a5_e*race+a6_e*disease+a7_e*insu+a8_e*x3)/sqrt(2))+1) - sqrt(2/pi)*a8_e*(r-1)*exp(-(a0_e+a1_e*income+a2_e*age+a3_e*gender+a4_e*educ+a5_e*race+a6_e*disease+a7_e*insu+a8_e*x3)^2/2)/(erf((a0_e+a1_e*income+a2_e*age+a3_e*gender+a4_e*educ+a5_e*race+a6_e*disease+a7_e*insu+a8_e*x3)/sqrt(2))-1) } #or not use erf but pnorm, same result.
      
      #two classic erros: 1. NA/inf error: log(1-pnorm) goes to -inf   2. trap, maximum updates: starting points are too few, need specify more
      cutoff <- (qnorm(1-10^-15)-a0_e-a1_e*income-a2_e*age-a3_e*gender-a4_e*educ-a5_e*race-a6_e*disease-a7_e*insu)/a8_e
      cutoff2 <- (qnorm(10^-15)-a0_e-a1_e*income-a2_e*age-a3_e*gender-a4_e*educ-a5_e*race-a6_e*disease-a7_e*insu)/a8_e
      lb <- -10   #original bound, could choose a relatively large number (to avoid trap error), since it will be further protected by cutoff to avoid NA/inf error.
      ub <- 10  #original bound
      if (a8_e<0) {
        range <- c(seq(max(ceiling(cutoff),lb),(max(ceiling(cutoff),lb)+min(floor(cutoff2),ub))/2, by=0.5), 
                   rev(seq(min(floor(cutoff2),ub), (max(ceiling(cutoff),lb)+min(floor(cutoff2),ub))/2, by=-0.5)))
        # range <- c(seq(max(cutoff,lb), 0, by=0.5), rev(seq(min(cutoff2,ub),0, by=-0.5)))
        mysample <- ars(n=v,f,fprima, m=length(range), x=range, lb=TRUE,xlb=cutoff,ub=TRUE,xub=cutoff2,
                        income=income, age=age, gende=gender, educ=educ, race=race, disease=disease, insu=insu,
                        b0_e=b0_e, b1_e=b1_e, b2_e=b2_e, b3_e=b3_e, b4_e=b4_e, b5_e=b5_e, b6_e=b6_e, sigma_e=sigma_e,
                        a0_e=a0_e, a1_e=a1_e, a2_e=a2_e, a3_e=a3_e, a4_e=a4_e, a5_e=a5_e, a6_e=a6_e, a7_e=a7_e, a8_e=a8_e, r=r)
      }else if (a8_e>0) {
        range <- c(seq(max(ceiling(cutoff2),lb),(max(ceiling(cutoff2),lb)+min(floor(cutoff),ub))/2, by=0.5), 
                   rev(seq(min(floor(cutoff),ub), (max(ceiling(cutoff2),lb)+min(floor(cutoff),ub))/2, by=-0.5)))
        # range <- c(seq(max(cutoff2,lb), 0, by=0.5), rev(seq(min(cutoff,ub),0, by=-0.5)))
        mysample <- ars(n=v,f,fprima, m=length(range), x=range, lb=TRUE,xlb=cutoff2,ub=TRUE,xub=cutoff,
                        income=income, age=age, gende=gender, educ=educ, race=race, disease=disease, insu=insu,
                        b0_e=b0_e, b1_e=b1_e, b2_e=b2_e, b3_e=b3_e, b4_e=b4_e, b5_e=b5_e, b6_e=b6_e, sigma_e=sigma_e,
                        a0_e=a0_e, a1_e=a1_e, a2_e=a2_e, a3_e=a3_e, a4_e=a4_e, a5_e=a5_e, a6_e=a6_e, a7_e=a7_e, a8_e=a8_e, r=r)
      }

      #newly sampled X3_r
      test_1[i,10:(10+v-1)] <- mysample
      remove(mysample)
    }
    list(test_1[i,])
    
  } #complete E-step 
  list_data <- Map(as.data.frame, res) 
  test_1  <- rbindlist(list_data, fill=TRUE)
  remove(list_data)
  remove(res)
  
  #M-step
  #re-format the data into long-form
  long1 <- subset(test_1, test_1$R==1, select=c(income, age, gender, educ, race, disease, insu, R,exp)) #without missing exp
  long1$w <- 1 
  colnames(long1)[9] <- c("exp_r")
  
  long2 <- data.frame(ID=1:(nrow(subset(test_1,test_1$R==0))*v),income=NA, age=NA, gender=NA, educ=NA, race=NA, disease=NA, insu=NA, R=0) #with missing X3
  temp <- subset(test_1, test_1$R==0, select=-c(exp,income, age, gender, educ, race, disease, insu, R))
  long2$exp_r <- as.vector(t(temp))
  for (i in 1:nrow(subset(test_1,test_1$R==0))){
    long2$educ[((i-1)*v+1):((i-1)*v+1+v-1)] <- subset(test_1, test_1$R==0, select=educ)$educ[i]
    long2$income[((i-1)*v+1):((i-1)*v+1+v-1)] <- subset(test_1, test_1$R==0, select=income)$income[i]
    long2$age[((i-1)*v+1):((i-1)*v+1+v-1)] <- subset(test_1, test_1$R==0, select=age)$age[i]
    long2$gender[((i-1)*v+1):((i-1)*v+1+v-1)] <- subset(test_1, test_1$R==0, select=gender)$gender[i]
    long2$race[((i-1)*v+1):((i-1)*v+1+v-1)] <- subset(test_1, test_1$R==0, select=race)$race[i]
    long2$disease[((i-1)*v+1):((i-1)*v+1+v-1)] <- subset(test_1, test_1$R==0, select=disease)$disease[i]
    long2$insu[((i-1)*v+1):((i-1)*v+1+v-1)] <- subset(test_1, test_1$R==0, select=insu)$insu[i]
  }
  long2 <- long2[,-1]
  long2$w <- 1/v
  
  long <- rbind(long1, long2) 
  
  
  #no need to update Xobs part under non-mixture
  #3) update outcome part: exp~age+gender+educ+race+disease+insu
  #note, for this update, we simply refit a weighted GLM reg 
  fit1 <- glm(formula=exp_r~age+gender+educ+race+disease+insu, data=long, weights=w)
  long$fitted3 <- fit1$fitted.values
  #update 
  para$b0_e[k+1] <- fit1[["coefficients"]][1]
  para$b1_e[k+1] <- fit1[["coefficients"]][2]
  para$b2_e[k+1] <- fit1[["coefficients"]][3]
  para$b3_e[k+1] <- fit1[["coefficients"]][4]
  para$b4_e[k+1] <- fit1[["coefficients"]][5]
  para$b5_e[k+1] <- fit1[["coefficients"]][6]
  para$b6_e[k+1] <- fit1[["coefficients"]][7]
  para$sigma_e[k+1] <- sqrt((sum((subset(long,long$w==1)$exp_r-subset(long,long$w==1)$fitted3)^2) + 
                             sum((1/v)*(subset(long,long$w==(1/v))$exp_r-subset(long,long$w==1/v)$fitted3)^2))/n)
  
  #4) update a0, a1, a2
  #note, for this update, we simply refit a weighted GLM reg 
  fit2 <- glm(formula=R~income+age+gender+educ+race+disease+insu+exp_r, long, family = binomial(link = "probit"), weights=w)
  
  para$a0_e[k+1]  <- fit2[["coefficients"]][1]
  para$a1_e[k+1]  <- fit2[["coefficients"]][2]
  para$a2_e[k+1]  <- fit2[["coefficients"]][3]
  para$a3_e[k+1]  <- fit2[["coefficients"]][4]
  para$a4_e[k+1]  <- fit2[["coefficients"]][5]
  para$a5_e[k+1]  <- fit2[["coefficients"]][6] 
  para$a6_e[k+1]  <- fit2[["coefficients"]][7]
  para$a7_e[k+1]  <- fit2[["coefficients"]][8]
  para$a8_e[k+1]  <- fit2[["coefficients"]][9]

  #convergence
  #squared distance for all key parameters are <0.5*10^-6 from previous iteration (equivalent to say most para has a 0.0002 change)
  if (k>100) {
    # if (sum((para[k,21:37]-para[k-1,21:37])^2) <10^(-7) |
    if (sum((para[k,21:37]-para[k-1,21:37])^2) <10^(-6) |
        (k==1000)) break
  }
  # kk <- rep(NA)
  # for (k in 2:500){
  #   kk[k] <- sum((para[k,21:37]-para[k-1,21:37])^2)
  # }
  # kk<10^(-7)
  # kk<2*10^(-5)
  
  
  remove(long)
  remove(long1)
  remove(long2)
  remove(temp)

  remove(fit1)
  remove(fit2)
  
} #one EM iteration completes. 















