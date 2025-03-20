library(ssmrob)
library(ars) #for ARS
library(pracma) #for the "erf" function
library(doParallel)
library (plyr)
library (dplyr)
library(data.table) 
library(sampleSelection)


registerDoParallel(cores=6)
data(MEPS2001)
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

para <- data.frame(iteration=1:2000,
                   prop=NA, #for binary cluster var "pik" (only consider two clusters)
                   #all below parameters are paired (for two clusters)
                   u0_a=NA, sd0_a=NA, u1_a=NA, sd1_a=NA, #age
                   u0_g=NA,  u1_g=NA, #gender
                   u0_r=NA,  u1_r=NA, #race
                   u0_e=NA, sd0_e=NA,  u1_e=NA, sd1_e=NA,  #edu
                   b00_inc=NA, b10_inc=NA, b20_inc=NA, b30_inc=NA, b40_inc=NA, sigma0_inc=NA,  #for income~age+gender+race+educ
                   b01_inc=NA, b11_inc=NA, b21_inc=NA, b31_inc=NA, b41_inc=NA, sigma1_inc=NA,  #for income~age+gender+race+educ
                   b00_ins=NA, b10_ins=NA, b01_ins=NA, b11_ins=NA, #for insu~income
                   b00_d=NA, b10_d=NA, b20_d=NA, b30_d=NA, b40_d=NA, b01_d=NA, b11_d=NA, b21_d=NA, b31_d=NA, b41_d=NA, #disease~age+gender+race+insu
                   b00_e=NA, b10_e=NA, b20_e=NA, b30_e=NA, b40_e=NA, b50_e=NA, b60_e=NA, sigma0_e=NA,  #for exp~age+gender+educ+race+disease+insu
                   b01_e=NA, b11_e=NA, b21_e=NA, b31_e=NA, b41_e=NA, b51_e=NA, b61_e=NA, sigma1_e=NA,  #for exp~age+gender+educ+race+disease+insu
                   a00_e=NA, a10_e=NA, a20_e=NA, a30_e=NA, a40_e=NA, a50_e=NA, a60_e=NA, a70_e=NA, a80_e=NA, 
                   a01_e=NA, a11_e=NA, a21_e=NA, a31_e=NA, a41_e=NA, a51_e=NA, a61_e=NA, a71_e=NA, a81_e=NA) #R~income+age+gender+educ+race+disease+insu

c1 <- coefficients(glm(formula=log(income)~age+gender+race+educ, data=test_1)) #note, income~logNormal
c2 <- coefficients(glm(formula=insu~income, data=test_1, family = binomial(link = "logit")))
c3 <- coefficients(glm(formula=disease~age+gender+race+insu, data=test_1,family = poisson))
c4 <- coefficients(glm(formula=exp~age+gender+educ+race+disease+insu, data=test_1))
c5 <- c(-0.100, 0.003, 0.113, 0.700, 0.066, -0.400, 0.860, 0.166, -0.116) #note, now we don't have full data to compute coeff, so we just use Miao's est

#for non-mixture, no iteration on below para except for exp~ and missing model (since no imputation for those var)
#i.e., only c4 and c5 need to be updated
#for mixture, all need to be updated
para[1,2] <- 0.88 #from Miao result
para[1,3] <- para[1,5] <- mean(test_1$age)
para[1,4] <- para[1,6] <- sd(test_1$age)
para[1,7] <- para[1,8] <- mean(test_1$gender)
para[1,9] <- para[1,10] <- mean(test_1$race)
para[1,11] <- para[1,13] <- mean(test_1$edu)
para[1,12] <- para[1,14] <- sd(test_1$edu)

para[1,15:19] <- c1
para[1,20] <- sqrt(var(lm(formula=log(income)~age+gender+race+educ, test_1)$residuals))
para[1,21:25] <- c1
para[1,26] <- sqrt(var(lm(formula=log(income)~age+gender+race+educ, test_1)$residuals)) + 0.5
para[1,27:28] <- para[1,29:30] <- c2
para[1,31:35] <- para[1,36:40] <- c3
para[1,41:47] <- para[1,49:55] <- c4 
para[1,47] <- 0.009  #from Miao result
para[1,55] <- -0.595 #from Miao result
para[1,48] <- 1.16  #from Miao result
para[1,56] <- 1.90  #from Miao result
para[1,57:65] <- para[1,66:74] <- c5
para[,74] <- 0  #MAR

n <- length(test_1$educ)
###################################################################################################################################
start.time <- Sys.time()

#EM starts here:
#E-step
test_t <- test_1
for (k in 1:1500) { #iteration for EM
  
  if (k %in% c(5, 10, 20, 30, 40) | k %% 50 ==0) {print(paste("EM-iter", k, "is running"))}
  test_1 <- test_t[,1:9]
  
  #E-step 
  #1. compute E(Z=z|Xobs,R,current theta) 2. compute E(logf(Xmis|Xobs) & E(logP(R|X) w.r.p. to Xmis|Xobs,R)
  #to obtain compute above two Expectations, especially when the E has no closed-form in general, 
  #we use (gibbs) sampling to obtain latent X3 & Z in order to compute Q function, per subject and per current EM iteration's parameter value
  # start.time <- Sys.time()
  res <- foreach(i=1:n) %dopar% {
    # for (i in 1:n) {
    
    #for those with missing exp, we need to sample both exp & Z
    if (test_1$R[i]==0){ #R=0 means missing
      
      #gibbs sampler for exp & Z (cluster) or f(exp,Z|Xobs,R;current parameter) per subject
      gibbs <- data.frame(iter=1:1000, exp=NA, Z=NA)
      gibbs[1,2:3] <- c(para$prop[k]*(para$b00_e[k]+para$b10_e[k]*test_1$age[i] +para$b20_e[k]*test_1$gender[i]+para$b30_e[k]*test_1$educ[i]+para$b40_e[k]*test_1$race[i]+para$b50_e[k]*test_1$disease[i]+para$b60_e[k]*test_1$insu[i])+
                        (1-para$prop[k])*(para$b01_e[k]+para$b11_e[k]*test_1$age[i] +para$b21_e[k]*test_1$gender[i]+para$b31_e[k]*test_1$educ[i]+para$b41_e[k]*test_1$race[i]+para$b51_e[k]*test_1$disease[i]+para$b61_e[k]*test_1$insu[i]), 1)    #initial values
      
      for (index in 1:1000) {
        #1) sample Z (we calculate a binary proportion for z=1, since we only consider two clusters), 
        #per P(Z|Xmis_t-1, Xobs, R; current para) = P(Z)*P(X1|Z)*P(X2|X1,Z)*...P(exp_t-1|.,Z)*P(R|X,Z) / sum of ...
        p1 <- para$prop[k] * #P(Z=1)
              dnorm(test_1$age[i], mean=para$u0_a[k], sd=para$sd0_a[k])* #f(age|Z=1)
              (para$u0_r[k]^test_1$race[i])*((1-para$u0_r[k])^(1-test_1$race[i]))* #f(race|Z=1)
              (para$u0_g[k]^test_1$gender[i])*((1-para$u0_g[k])^(1-test_1$gender[i]))* #f(gender|Z=1)
              dnorm(test_1$educ[i], mean=para$u0_e[k], sd=para$sd0_e[k]) * #f(educ|Z=1)
              dlnorm(test_1$income[i], mean=para$b00_inc[k]+para$b10_inc[k]*test_1$age[i] +para$b20_inc[k]*test_1$gender[i]+para$b30_inc[k]*test_1$race[i]+para$b40_inc[k]*test_1$educ[i],
                                       sd=para$sigma0_inc[k])*  #f(income|.,Z=1) -lognormal
              (1/(1+exp(-para$b00_ins[k]-para$b10_ins[k]*test_1$income[i]))^test_1$insu[i])*((1-1/(1+exp(-para$b00_ins[k]-para$b10_ins[k]*test_1$income[i])))^(1-test_1$insu[i]))*  #f(insurance|.,Z=1)
              dpois(test_1$disease[i], lambda=exp(para$b00_d[k]+para$b10_d[k]*test_1$age[i] +para$b20_d[k]*test_1$gender[i]+para$b30_d[k]*test_1$race[i]+para$b40_d[k]*test_1$insu[i]))*  #f(disease|.,Z=1)
              dnorm(gibbs[index,2], mean=para$b00_e[k]+para$b10_e[k]*test_1$age[i] +para$b20_e[k]*test_1$gender[i]+para$b30_e[k]*test_1$educ[i]+para$b40_e[k]*test_1$race[i]+para$b50_e[k]*test_1$disease[i]+para$b60_e[k]*test_1$insu[i], 
                                    sd=para$sigma0_e[k]) * #f(exp|.,Z=1)
              (1-pnorm(para$a00_e[k]+para$a10_e[k]*test_1$income[i] +para$a20_e[k]*test_1$age[i]+para$a30_e[k]*test_1$gender[i]+para$a40_e[k]*test_1$educ[i]+para$a50_e[k]*test_1$race[i]+para$a60_e[k]*test_1$disease[i]+para$a70_e[k]*test_1$insu[i]+para$a80_e[k]*gibbs[index,2]))  #P(R=0|.,Z=1)
        p2 <-  (1-para$prop[k]) * #P(Z=0)
               dnorm(test_1$age[i], mean=para$u1_a[k], sd=para$sd1_a[k])* #f(age|Z=0)
               (para$u1_r[k]^test_1$race[i])*((1-para$u1_r[k])^(1-test_1$race[i]))* #f(race|Z=0)
               (para$u1_g[k]^test_1$gender[i])*((1-para$u1_g[k])^(1-test_1$gender[i]))* #f(gender|Z=0)
               dnorm(test_1$educ[i], mean=para$u1_e[k], sd=para$sd1_e[k]) * #f(educ|Z=0)
               dlnorm(test_1$income[i], mean=para$b01_inc[k]+para$b11_inc[k]*test_1$age[i] +para$b21_inc[k]*test_1$gender[i]+para$b31_inc[k]*test_1$race[i]+para$b41_inc[k]*test_1$educ[i],
                                        sd=para$sigma1_inc[k])*  #f(income|.,Z=0) 
               (1/(1+exp(-para$b01_ins[k]-para$b11_ins[k]*test_1$income[i]))^test_1$insu[i])*((1-1/(1+exp(-para$b01_ins[k]-para$b11_ins[k]*test_1$income[i])))^(1-test_1$insu[i]))*  #f(insurance|.,Z=0)
               dpois(test_1$disease[i], lambda=exp(para$b00_d[k]+para$b10_d[k]*test_1$age[i] +para$b20_d[k]*test_1$gender[i]+para$b30_d[k]*test_1$race[i]+para$b40_d[k]*test_1$insu[i]))*  #f(disease|.,Z=0)
               dnorm(gibbs[index,2], mean=para$b01_e[k]+para$b11_e[k]*test_1$age[i] +para$b21_e[k]*test_1$gender[i]+para$b31_e[k]*test_1$educ[i]+para$b41_e[k]*test_1$race[i]+para$b51_e[k]*test_1$disease[i]+para$b61_e[k]*test_1$insu[i], 
                                     sd=para$sigma1_e[k]) * #f(exp|.,Z=0)
               (1-pnorm(para$a01_e[k]+para$a11_e[k]*test_1$income[i] +para$a21_e[k]*test_1$age[i]+para$a31_e[k]*test_1$gender[i]+para$a41_e[k]*test_1$educ[i]+para$a51_e[k]*test_1$race[i]+para$a61_e[k]*test_1$disease[i]+para$a71_e[k]*test_1$insu[i]+para$a81_e[k]*gibbs[index,2]))  #P(R=0|.,Z=0)
        
        # p2 <- (1-para$prop[k]) *   
        #   dnorm(test_1$X1[i], mean=para$u01[k], sd=para$sd01[k]) * 
        #   dnorm(test_1$X2[i], mean=para[k,10]+para[k,11]*test_1$X1[i], sd=para$sd1[k]) *
        #   dnorm(gibbs[index,2], mean=para[k,16]+para[k,17]*test_1$X2[i], sd=para$sigma1[k]) *
        #   (1-pnorm(para[k,22]+0*gibbs[index,2]+para[k,24]*test_1$X2[i]))    #for MAR
        prop1 <- p1/(p1+p2) #should be ~1 (1: labeled into 1st (true) cluster; 0: 2nd (wrong) cluster)
        
        #newly sampled Z_r for gibbs, then sample exp
        gibbs[index+1, 3] <- rbinom(1,1,prop1)
        
        #2) sample exp, we use ADS (since P(exp|age+gender+educ+race+disease+insu,R) prop to 
        #                                 P(exp|age+gender+educ+race+disease+insu  )*P(R|exp,income+age+gender+educ+race+disease+insu))
        #note, the ADS here is conditioned on previously sampled Z_r (1 or 0), hence two separate ADS based on Z.
        if (gibbs[index+1, 3] == 1){
          #current EM iteration parameter value for cluster 1
          b0_e <- para[k,41]
          b1_e <- para[k,42]
          b2_e <- para[k,43]
          b3_e <- para[k,44]
          b4_e <- para[k,45]
          b5_e <- para[k,46]
          b6_e <- para[k,47]
          sigma_e <- para[k,48]
          a0_e  <- para[k,57]
          a1_e  <- para[k,58] 
          a2_e  <- para[k,59] 
          a3_e  <- para[k,60] 
          a4_e  <- para[k,61] 
          a5_e  <- para[k,62] 
          a6_e  <- para[k,63] 
          a7_e  <- para[k,64]
          a8_e  <- para[k,65]
          
          #draw "1" sample for missing X3
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
            mysample <- ars(n=1,f,fprima, m=length(range), x=range, lb=TRUE,xlb=cutoff,ub=TRUE,xub=cutoff2,
                            income=income, age=age, gende=gender, educ=educ, race=race, disease=disease, insu=insu,
                            b0_e=b0_e, b1_e=b1_e, b2_e=b2_e, b3_e=b3_e, b4_e=b4_e, b5_e=b5_e, b6_e=b6_e, sigma_e=sigma_e,
                            a0_e=a0_e, a1_e=a1_e, a2_e=a2_e, a3_e=a3_e, a4_e=a4_e, a5_e=a5_e, a6_e=a6_e, a7_e=a7_e, a8_e=a8_e, r=r)
          }else if (a8_e>0) {
            range <- c(seq(max(ceiling(cutoff2),lb),(max(ceiling(cutoff2),lb)+min(floor(cutoff),ub))/2, by=0.5), 
                       rev(seq(min(floor(cutoff),ub), (max(ceiling(cutoff2),lb)+min(floor(cutoff),ub))/2, by=-0.5)))
            # range <- c(seq(max(cutoff2,lb), 0, by=0.5), rev(seq(min(cutoff,ub),0, by=-0.5)))
            mysample <- ars(n=1,f,fprima, m=length(range), x=range, lb=TRUE,xlb=cutoff2,ub=TRUE,xub=cutoff,
                            income=income, age=age, gende=gender, educ=educ, race=race, disease=disease, insu=insu,
                            b0_e=b0_e, b1_e=b1_e, b2_e=b2_e, b3_e=b3_e, b4_e=b4_e, b5_e=b5_e, b6_e=b6_e, sigma_e=sigma_e,
                            a0_e=a0_e, a1_e=a1_e, a2_e=a2_e, a3_e=a3_e, a4_e=a4_e, a5_e=a5_e, a6_e=a6_e, a7_e=a7_e, a8_e=a8_e, r=r)
          }
          #newly sampled X3_r
          gibbs[index+1, 2] <- mysample
          remove(mysample)
          
        }else if (gibbs[index+1, 3] == 0) {
          # for MAR: note, in this case, P(X3|X1,X2,R,Z) just prop to P(X3|X2,Z), no ARS needed.
          # i.e., P(exp|age+gender+educ+race+disease+insu,R) prop to P(exp|age+gender+educ+race+disease+insu  )
          b0_e <- para[k,49]
          b1_e <- para[k,50]
          b2_e <- para[k,51]
          b3_e <- para[k,52]
          b4_e <- para[k,53]
          b5_e <- para[k,54]
          b6_e <- para[k,55]
          sigma_e <- para[k,56]
          
          #draw "1" sample for missing X3
          income <- test_1$income[i]
          age <- test_1$age[i]
          gender <- test_1$gender[i]
          educ <- test_1$educ[i]
          race <- test_1$race[i]
          disease <- test_1$disease[i]
          insu <- test_1$insu[i]

          gibbs[index+1, 2] <- rnorm(1, mean=b0_e+b1_e*age+b2_e*gender+b3_e*educ+b4_e*race+b5_e*disease+b6_e*insu, sd=sigma_e)
        }
      }
      
      #obtain (Z,X3)|X1,X2,R, current para, for E-step
      #note, to compute the parenthesis part expectation (w.r.p to X3|Z,Xobs), we will use two different subsamples (per Z) from above marginal sample of X3
      #i.e., for M-step, each subject will be used twice in two clusters, each with a double-weight (=marginal pi_k*weight from subsample) 
      ss <- seq(201,1000,by=2)  #burn-in=1000, thinning=2
      
      test_1[i,10] <- sum(gibbs$Z[ss] == 1)/length(ss)  #for P(Z|Xobs,R;current para), 1st weight
      test_1[i,11] <- sum(gibbs$Z[ss] == 1)
      colnames(test_1)[10] <- c("prop1")
      colnames(test_1)[11] <- c("n1")
      
      sub1 <- subset(gibbs$exp, gibbs$Z==1 & gibbs$iter %in% ss) #for P(Xmis|Xobs,R,Z;current para), its sample size =1/2nd weight
      sub2 <- subset(gibbs$exp, gibbs$Z==0 & gibbs$iter %in% ss) #for P(Xmis|Xobs,R,Z;current para), its sample size =1/2nd weight
      if (length(sub1)>0 & length(sub2)>0){
        test_1[i,12:(10+length(sub1)+1)] <- sub1
        test_1[i,(10+length(sub1)+2):(10+length(sub1)+1+length(sub2))] <- sub2
      }else if (length(sub1)==0 & length(sub2)>0){
        test_1[i,(10+length(sub1)+2):(10+length(sub1)+1+length(sub2))] <- sub2
      }else if (length(sub2)==0 & length(sub1)>0){
        test_1[i,12:(10+length(sub1)+1)] <- sub1
      }
      remove(sub1)
      remove(sub2)
      remove(gibbs)
      
    }else if (test_1$R[i]==1){ #for those without missing X3, we only need to consider Z (and no gibbs)
      
      #calculate exact P(Z|X,R;current para) (only difference with above gibbs is that now Xmis is fixed, no gibbs on Xmis)
      #P(Z|Xmis, Xobs, R; current para) = P(Z)*P(X1|Z)*P(X2|X1,Z)*P(X3|X1,X2,Z)*P(R|X) / sum of ...
      #note, for this P(Z|X,R;current para), when marginal P(Z=1) is ~1, P(Z=1|X,R;current para) will almost surely dominate by term P(Z=1), 
      #i.e., almost all data (even for Z=0) will be labeled into cluster 1 just because of this term P(Z=1) >> P(Z=0), 
      #this phenomnon will be attenuated only when two clusters have extremely different P(X|Z) and/or P(R|X,Z) to reverse this term's effect.
      #compute P(X3,X1,X2,Z=z,R;current para) for Z=1,0; then divided by sum or P(X3,X1,X2,R;current para).
      p1 <- para$prop[k] * #P(Z=1)
            dnorm(test_1$age[i], mean=para$u0_a[k], sd=para$sd0_a[k])* #f(age|Z=1)
            (para$u0_r[k]^test_1$race[i])*((1-para$u0_r[k])^(1-test_1$race[i]))* #f(race|Z=1)
            (para$u0_g[k]^test_1$gender[i])*((1-para$u0_g[k])^(1-test_1$gender[i]))* #f(gender|Z=1)
            dnorm(test_1$educ[i], mean=para$u0_e[k], sd=para$sd0_e[k]) * #f(educ|Z=1)
            dlnorm(test_1$income[i], mean=para$b00_inc[k]+para$b10_inc[k]*test_1$age[i] +para$b20_inc[k]*test_1$gender[i]+para$b30_inc[k]*test_1$race[i]+para$b40_inc[k]*test_1$educ[i],
                                     sd=para$sigma0_inc[k])*  #f(income|.,Z=1) -lognormal
            (1/(1+exp(-para$b00_ins[k]-para$b10_ins[k]*test_1$income[i]))^test_1$insu[i])*((1-1/(1+exp(-para$b00_ins[k]-para$b10_ins[k]*test_1$income[i])))^(1-test_1$insu[i]))*  #f(insurance|.,Z=1)
            dpois(test_1$disease[i], lambda=exp(para$b00_d[k]+para$b10_d[k]*test_1$age[i] +para$b20_d[k]*test_1$gender[i]+para$b30_d[k]*test_1$race[i]+para$b40_d[k]*test_1$insu[i]))*  #f(disease|.,Z=1)
            dnorm(test_1$exp[i], mean=para$b00_e[k]+para$b10_e[k]*test_1$age[i] +para$b20_e[k]*test_1$gender[i]+para$b30_e[k]*test_1$educ[i]+para$b40_e[k]*test_1$race[i]+para$b50_e[k]*test_1$disease[i]+para$b60_e[k]*test_1$insu[i], 
                                 sd=para$sigma0_e[k]) * #f(exp|.,Z=1)
            (pnorm(para$a00_e[k]+para$a10_e[k]*test_1$income[i] +para$a20_e[k]*test_1$age[i]+para$a30_e[k]*test_1$gender[i]+para$a40_e[k]*test_1$educ[i]+para$a50_e[k]*test_1$race[i]+para$a60_e[k]*test_1$disease[i]+para$a70_e[k]*test_1$insu[i]+para$a80_e[k]*test_1$exp[i]))  #P(R=1|.,Z=1)
      p2 <-  (1-para$prop[k]) * #P(Z=0)
             dnorm(test_1$age[i], mean=para$u1_a[k], sd=para$sd1_a[k])* #f(age|Z=0)
             (para$u1_r[k]^test_1$race[i])*((1-para$u1_r[k])^(1-test_1$race[i]))* #f(race|Z=0)
             (para$u1_g[k]^test_1$gender[i])*((1-para$u1_g[k])^(1-test_1$gender[i]))* #f(gender|Z=0)
             dnorm(test_1$educ[i], mean=para$u1_e[k], sd=para$sd1_e[k]) * #f(educ|Z=0)
             dlnorm(test_1$income[i], mean=para$b01_inc[k]+para$b11_inc[k]*test_1$age[i] +para$b21_inc[k]*test_1$gender[i]+para$b31_inc[k]*test_1$race[i]+para$b41_inc[k]*test_1$educ[i],
                                      sd=para$sigma1_inc[k])*  #f(income|.,Z=0) 
             (1/(1+exp(-para$b01_ins[k]-para$b11_ins[k]*test_1$income[i]))^test_1$insu[i])*((1-1/(1+exp(-para$b01_ins[k]-para$b11_ins[k]*test_1$income[i])))^(1-test_1$insu[i]))*  #f(insurance|.,Z=0)
             dpois(test_1$disease[i], lambda=exp(para$b00_d[k]+para$b10_d[k]*test_1$age[i] +para$b20_d[k]*test_1$gender[i]+para$b30_d[k]*test_1$race[i]+para$b40_d[k]*test_1$insu[i]))*  #f(disease|.,Z=0)
             dnorm(test_1$exp[i], mean=para$b01_e[k]+para$b11_e[k]*test_1$age[i] +para$b21_e[k]*test_1$gender[i]+para$b31_e[k]*test_1$educ[i]+para$b41_e[k]*test_1$race[i]+para$b51_e[k]*test_1$disease[i]+para$b61_e[k]*test_1$insu[i], 
                                  sd=para$sigma1_e[k]) * #f(exp|.,Z=0)
             (pnorm(para$a01_e[k]+para$a11_e[k]*test_1$income[i] +para$a21_e[k]*test_1$age[i]+para$a31_e[k]*test_1$gender[i]+para$a41_e[k]*test_1$educ[i]+para$a51_e[k]*test_1$race[i]+para$a61_e[k]*test_1$disease[i]+para$a71_e[k]*test_1$insu[i]+para$a81_e[k]*test_1$exp[i]))  #P(R=1|.,Z=0)
      # p2 <- (1-para$prop[k]) * 
      #   dnorm(test_1$X1[i], mean=para$u01[k], sd=para$sd01[k]) * 
      #   dnorm(test_1$X2[i], mean=para[k,10]+para[k,11]*test_1$X1[i], sd=para$sd1[k]) *
      #   dnorm(test_1$X3_r[i], mean=para[k,16]+para[k,17]*test_1$X2[i], sd=para$sigma1[k]) *
      #   pnorm(para[k,22]+0*test_1$X3_r[i]+para[k,24]*test_1$X2[i])   #for MAR
      
      #obtain P(Z=1|Xmis,Xobs,R;current para), 1st (and sole) weight 
      test_1[i,10] <-  p1/(p1+p2) 
      colnames(test_1)[10] <- c("prop1")
      
    }
    list(test_1[i,])
    
  } #complete E-step 
  list_data <- Map(as.data.frame, res) 
  test_1  <- rbindlist(list_data, fill=TRUE)
  remove(list_data)
  remove(res)
  
  
  #M-step
  #parameters to be updated: 
  #1) pi_K  
  #2) u & sd for X1  
  #3) b0 b1 & sd for X2|X1   
  #4) beta0 beta1 & sigma for X3|X2
  #5) a0 a1 a2 for R|X2,X3
 
  #re-format the data into long-form
  test_1$ID <- seq(1,3323,by=1)
  long11 <- long12 <- subset(test_1, test_1$R==1, select=c(ID, income, age, gender, educ, race, disease, insu, R, prop1, exp)) #CC
  long11$cluster <- 1  #component with bigger N
  long12$prop1 <- 1- long12$prop1
  long12$cluster <- 0  #component with smaller N
  long1 <- rbind(long11, long12) 
  long1 <- long1[order(long1$ID),]
  long1$w <- 1 
  long1 <- long1[,-1]
  remove(long11)
  remove(long12)
  colnames(long1)[10] <- c("exp_r")

  long2 <- data.frame(ID=1:(nrow(subset(test_1,test_1$R==0))*400),income=NA, age=NA, gender=NA, educ=NA, race=NA, disease=NA, insu=NA, R=0, prop1=NA, exp_r=NA, cluster=NA, w=NA) #with missing X3
  temp <- subset(test_1, test_1$R==0, select=-c(ID,exp,income, age, gender, educ, race, disease, insu, R, prop1, n1))
  long2$exp_r <- as.vector(t(temp))
  for (i in 1:nrow(subset(test_1,test_1$R==0))){
    long2$educ[((i-1)*400+1):((i-1)*400+1+400-1)] <- subset(test_1, test_1$R==0, select=educ)$educ[i]
    long2$income[((i-1)*400+1):((i-1)*400+1+400-1)] <- subset(test_1, test_1$R==0, select=income)$income[i]
    long2$age[((i-1)*400+1):((i-1)*400+1+400-1)] <- subset(test_1, test_1$R==0, select=age)$age[i]
    long2$gender[((i-1)*400+1):((i-1)*400+1+400-1)] <- subset(test_1, test_1$R==0, select=gender)$gender[i]
    long2$race[((i-1)*400+1):((i-1)*400+1+400-1)] <- subset(test_1, test_1$R==0, select=race)$race[i]
    long2$disease[((i-1)*400+1):((i-1)*400+1+400-1)] <- subset(test_1, test_1$R==0, select=disease)$disease[i]
    long2$insu[((i-1)*400+1):((i-1)*400+1+400-1)] <- subset(test_1, test_1$R==0, select=insu)$insu[i]

    if (subset(test_1,test_1$R==0)$n1[i] ==0){  #all sampled Z are in cluster 0
      long2$prop1[(((i-1)*400+1)):((i-1)*400+400)]  <- 1
      long2$w[(((i-1)*400+1)):((i-1)*400+400)]  <- 1/400
      long2$cluster[(((i-1)*400+1)):((i-1)*400+400)]  <- 0
    }else if (subset(test_1,test_1$R==0)$n1[i] ==400){ #all sampled Z are in cluster 1
      long2$prop1[((i-1)*400+1):(((i-1)*400+1)+400-1)]  <-  1
      long2$w[((i-1)*400+1):(((i-1)*400+1)+400-1)]  <-  1/400
      long2$cluster[((i-1)*400+1):(((i-1)*400+1)+400-1)]  <-  1
    }else {
      long2$prop1[((i-1)*400+1):(((i-1)*400+1)+subset(test_1,test_1$R==0)$n1[i]-1)]  <-  subset(test_1,test_1$R==0)$prop1[i]
      long2$prop1[(((i-1)*400+1)+subset(test_1,test_1$R==0)$n1[i]):((i-1)*400+400)]  <- 1-subset(test_1,test_1$R==0)$prop1[i]
      long2$w[((i-1)*400+1):(((i-1)*400+1)+subset(test_1,test_1$R==0)$n1[i]-1)]  <-  1/subset(test_1,test_1$R==0)$n1[i]
      long2$w[(((i-1)*400+1)+subset(test_1,test_1$R==0)$n1[i]):((i-1)*400+400)]  <- 1/(400-subset(test_1,test_1$R==0)$n1[i])
      long2$cluster[((i-1)*400+1):(((i-1)*400+1)+subset(test_1,test_1$R==0)$n1[i]-1)]  <-  1
      long2$cluster[(((i-1)*400+1)+subset(test_1,test_1$R==0)$n1[i]):((i-1)*400+400)]  <-  0
    }
  }
  long2 <- long2[,-1]
  remove(temp)
  
  long <- rbind(long1, long2) 
  long$w2 <- long$w*long$prop1
  for (i in 1:length(long$exp_r)){
    if (is.na(long$w2[i])=='TRUE') {long$w2[i] <- 0}
  }
  # sum(long$w2) =N
  
  #1) update marginal pi_k
  #note, for pi_k, only the 1st weight is needed (E(I(Zi=k)) w.r.p to Z|Xobs,R), and it has a closed-form MLE solution.
  #for those without missing, the weight is exact; for those with missing, the weight is computed via gibbs sampler.
  para$prop[k+1] <- sum(test_1$prop1)/n   #n=(sum(test_1$prop1)+sum(1-test_1$prop1))
  # summary(lm(formula=prop1~1, data=test_1))
  
  #2) update X
  #note, here we will use two weights (one from E(I(Zi=k)), one from gibbs sampler for X3|X2,current para)
  #2.1 update 4 baseline X para 
  d1 <- subset(long,long$cluster==1)
  d2 <- subset(long,long$cluster==0)
  fit_age1 <- lm(formula=age~1, data=d1, weights=w2)
  fit_age2 <- lm(formula=age~1, data=d2, weights=w2)
  d1$fittedage1 <- fit_age1$fitted.values
  d2$fittedage1 <- fit_age2$fitted.values
  
  para[k+1,3] <- fit_age1[["coefficients"]]
  para[k+1,5] <- fit_age2[["coefficients"]]
  para[k+1,4] <- sqrt(sum((d1$age-d1$fittedage1)^2*d1$w2) /sum(d1$w2))
  para[k+1,6] <- sqrt(sum((d2$age-d2$fittedage1)^2*d2$w2) /sum(d2$w2))
  
  fit_educ1 <- lm(formula=educ~1, data=d1, weights=w2)
  fit_educ2 <- lm(formula=educ~1, data=d2, weights=w2)
  d1$fittededuc1 <- fit_educ1$fitted.values
  d2$fittededuc1 <- fit_educ2$fitted.values
  para[k+1,11] <- fit_educ1[["coefficients"]]
  para[k+1,13] <- fit_educ2[["coefficients"]]
  para[k+1,12] <- sqrt(sum((d1$educ-d1$fittededuc1)^2*d1$w2) /sum(d1$w2))
  para[k+1,14] <- sqrt(sum((d2$educ-d2$fittededuc1)^2*d2$w2) /sum(d2$w2))
  
  para[k+1,7] <- sum(d1$gender*d1$w2)/sum(d1$w2)
  para[k+1,8] <- sum(d2$gender*d2$w2)/sum(d2$w2)
  
  para[k+1,9] <- sum(d1$race*d1$w2)/sum(d1$w2)
  para[k+1,10] <- sum(d2$race*d2$w2)/sum(d2$w2)
  
  para[k+1,11] <- para[1,13] <- mean(test_1$edu)
  para[k+1,12] <- para[1,14] <- sd(test_1$edu)
  
  #3) update rest of 4 X (conditional reg)
  fit_inc1 <- lm(formula=log(income)~age+gender+race+educ, data=d1, weights=w2)
  fit_inc2 <- lm(formula=log(income)~age+gender+race+educ, data=d2, weights=w2)
  d1$fittedinc1 <- fit_inc1$fitted.values
  d2$fittedinc1 <- fit_inc2$fitted.values
  para[k+1,15:19] <- fit_inc1[["coefficients"]]
  para[k+1,20] <- sqrt(sum((log(d1$income)-(d1$fittedinc1))^2*d1$w2) /sum(d1$w2))
  para[k+1,21:25] <- fit_inc2[["coefficients"]]
  para[k+1,26] <- sqrt(sum((log(d2$income)-(d2$fittedinc1))^2*d2$w2) /sum(d2$w2))
  
  fit_insu1 <- glm(formula=insu~income, data=d1, family = binomial(link = "logit"), weights=w2)
  fit_insu2 <- glm(formula=insu~income, data=d2, family = binomial(link = "logit"), weights=w2)
  para[k+1,27:28] <- fit_insu1[["coefficients"]]
  para[k+1,29:30] <- fit_insu2[["coefficients"]]
  
  fit_disease1 <- glm(formula=disease~age+gender+race+insu, family = poisson, data=d1, weights=w2)
  fit_disease2 <- glm(formula=disease~age+gender+race+insu, family = poisson, data=d2, weights=w2)
  para[k+1,31:35] <- fit_disease1[["coefficients"]]
  para[k+1,36:40] <- fit_disease2[["coefficients"]]
  
  fit_exp1 <- lm(formula=exp_r~age+gender+educ+race+disease+insu, data=d1, weights=w2)
  fit_exp2 <- lm(formula=exp_r~age+gender+educ+race+disease+insu, data=d2, weights=w2)
  d1$fittedexp1 <- fit_exp1$fitted.values
  d2$fittedexp1 <- fit_exp2$fitted.values
  para[k+1,41:47] <- fit_exp1[["coefficients"]]
  para[k+1,49:55] <- fit_exp2[["coefficients"]]
  para[k+1,48] <- sqrt(sum((d1$exp_r-d1$fittedexp1)^2*d1$w2) /sum(d1$w2))
  para[k+1,56] <- sqrt(sum((d2$exp_r-d2$fittedexp1)^2*d2$w2) /sum(d2$w2))
  
  #4) update missing model R|X
  fit_R1 <- glm(formula=R~income+age+gender+educ+race+disease+insu+exp_r, d1, family = binomial(link = "probit"), weights=w2)
  fit_R2 <- glm(formula=R~income+age+gender+educ+race+disease+insu , d2, family = binomial(link = "probit"), weights=w2)
  para[k+1,57:65] <- fit_R1[["coefficients"]]
  para[k+1,66:73] <- fit_R2[["coefficients"]]
  
  remove(d1)
  remove(d2)
  remove(long)
  remove(long1)
  remove(long2)
  
  remove(fit_age1)
  remove(fit_age2)
  remove(fit_educ1)
  remove(fit_educ2)
  remove(fit_inc1)
  remove(fit_inc2)
  remove(fit_insu1)
  remove(fit_insu2)
  remove(fit_disease1)
  remove(fit_disease2)
  remove(fit_exp1)
  remove(fit_exp2)
  remove(fit_R1)
  remove(fit_R2)
  
  #convergence
  #squared distance for all parameters are <10^-5 from previous iteration (equivalent to say the key 17 para has a 0.0002 change)
  if (k>200) {
    if (sum((para[k,2:74]-para[k-1,2:74])^2) <10^(-6)|
        sum((para[k,-1]-para[k-1,-1])^2) <10^(-6) |
        (k==2000))  break
  }

  
} #one EM iteration completes. 
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

####################################################################################################################################


saveRDS(test_1, file="/Users/riceball/Desktop/PhD/Paper/Paper 1 - Mixture MNAR/realcase_2 component (MNAR+MAR)_test_1.csv")
saveRDS(para, file="/Users/riceball/Desktop/PhD/Paper/Paper 1 - Mixture MNAR/realcase_2 component_(MNAR+MAR)_para.csv")
saveRDS(d1, file="/Users/riceball/Desktop/PhD/Paper/Paper 1 - Mixture MNAR/realcase_2 component_(MNAR+MAR)_d1.csv")
saveRDS(d2, file="/Users/riceball/Desktop/PhD/Paper/Paper 1 - Mixture MNAR/realcase_2 component_(MNAR+MAR)_d2.csv")

#pateint characteristics
for (i in 1:3323){
  if (test_1$prop1[i] >0.5) {
    test_1$cluster[i] <- 0
    }else { test_1$cluster[i] <- 1}
}

table(test_1$cluster)
# install.packages("purrr")
library("purrr")
descrip <- test_1[,c(1:9,413)]
descrip %>% split(.$cluster) %>%map(summary)
sd(subset(test_1, cluster==0)$age)
sd(subset(test_1, cluster==1)$age)







