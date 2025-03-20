library(ars) #for ARS
library(pracma) #for the "erf" function
library(doParallel)
library (plyr)
library(data.table) 

set.seed(12345)

paraa <- list()
testt <- list() 

# saveRDS(paraa, file="/Users/riceball/Desktop/PhD/Paper/Paper 1 - Mixture MNAR/paraa_mixture+mixture(10 miss)_n=500.csv")
# saveRDS(testt, file="/Users/riceball/Desktop/PhD/Paper/Paper 1 - Mixture MNAR/testt_mixture+mixture(10 miss)_n=500.csv")
paraa <- readRDS("/Users/riceball/Desktop/PhD/Paper/Paper 1 - Mixture MNAR/paraa_mixture+mixture(10 miss)_n=500.csv")
testt <- readRDS("/Users/riceball/Desktop/PhD/Paper/Paper 1 - Mixture MNAR/testt_mixture+mixture(10 miss)_n=500.csv")


saveRDS(paraa, file="/Users/riceball/Desktop/PhD/Paper/Paper 1 - Mixture MNAR/paraa_mixture+mixture(30 miss)_n=1000.csv")
saveRDS(testt, file="/Users/riceball/Desktop/PhD/Paper/Paper 1 - Mixture MNAR/testt_mixture+mixture(30 miss)_n=1000.csv")
paraa <- readRDS("/Users/riceball/Desktop/PhD/Paper/Paper 1 - Mixture MNAR/paraa_mixture+mixture(30 miss)_n=1000.csv")
testt <- readRDS("/Users/riceball/Desktop/PhD/Paper/Paper 1 - Mixture MNAR/testt_mixture+mixture(30 miss)_n=1000.csv")


# parallel::detectCores() # number of cores on this computer
registerDoParallel(cores=4)

###################################################################################################################################
#note, I tested different combination of N*prop, even when N=2000, if prop=90%, the cluster 2 estimation is still biased.
n <- 1000 #size of each simulated dataset
prop <- 0.9 #true mixing weight


start.time <- Sys.time()
for (t in 1:1000){ #generate 100 samples
  if (t %% 1 ==0) {print(paste("sample", t, "is being processed"))}
  
  #note, if the true is a non-mixture, no doubt the estimate for cluster 1 is "biased", this is a trade-off of bias vs. uncertainty of missingness.
  #true data 1: no mixture
  # test_1 <- data.frame(ID=1:n, X3=NA, X3_r=NA, X2=NA, X1=NA, R=NA, cluster=NA)
  # test_1$X1 <- rnorm(n,mean=0,sd=0.5)
  # test_1$X2 <- rnorm(n,mean=0.8+0.5*test_1$X1,sd=0.5) #note, here rather directly specifying X2, we made it into a function of X1 so that afterwards when using P(X1,X2|Z) we could break it into two parts.
  # test_1$X3 <- rnorm(n,mean=0.8+0.8*test_1$X2,sd=sqrt(1)) 
  # test_1$R <- rbinom(n,1, pnorm(1 + 0.8*test_1$X3 - 0.5*test_1$X2))
  # for (i in 1:n){
  #   if (test_1$R[i] == 1) {test_1$X3_r[i] <- test_1$X3[i]}
  # }
  # table(test_1$R)
  # summary(test_1$X1)
  # sd(test_1$X1)
  # lm(formula=X2~X1, test_1)
  # lm(formula=X3~X2, test_1)
  # glm(formula=R~X3+X2, test_1, family = binomial(link = "probit"))
  
  #true data 2: mixture
  test_1 <- data.frame(ID=1:n, X3=NA, X3_r=NA, X2=NA, X1=NA, R=NA, cluster=NA)
  test_1$cluster <- rbinom(n,1,prop)
  for (i in 1:n){
    if (test_1$cluster[i] == 1) { #cluster 1 ~true  (order: 1 >> 0)
      test_1$X1[i] <- rnorm(1,mean=0,sd=0.5) #marginal X1 has a ordering
      test_1$X2[i] <- rnorm(1,mean=0.8+0.5*test_1$X1[i],sd=0.5) #b2 must <= b2' (coef for X1, i.e., 0.5<1); sd2 <= sd2'(0.5<=0.5); also, marginal X2 needs a ordering
      test_1$X3[i] <- rnorm(1,mean=0.8+0.8*test_1$X2[i],sd=sqrt(1))  #a2 must <= a2' (coef for X2, i.e., 0.8<=0.8); sd3 <= sd3'(1<=1);also, marginal X3 needs a ordering

      # test_1$R3_s[i] <- 1.0 + 0.8*test_1$X3[i] - 0.5*test_1$X2[i] +rnorm(1,mean=0,sd=1) #coef for X3 needs to be different
      # # test_1$R[i] <- ifelse(test_1$R3_s[i]>0,1,0) #R=1 means observed.
      # test_1$R[i] <- rbinom(1,1, pnorm(1 + 0.8*test_1$X3[i] - 0.5*test_1$X2[i]))
      test_1$R[i] <- rbinom(1,1, pnorm(1 + 0.8*test_1$X3[i] - 0.5*test_1$X2[i])) #10% miss
      # test_1$R[i] <- rbinom(1,1, pnorm(2.5 - 0.9*test_1$X3[i] - 0.5*test_1$X2[i] ))  #30% miss
    } else { #cluster 2 ～wrong
      test_1$X1[i] <- rnorm(1,mean=1.5,sd=0.5)
      test_1$X2[i] <- rnorm(1,mean=2.0+0.5*test_1$X1[i],sd=0.5)  #to satisfy condition 1, when sigma for X1 is same, we let beta for X1 being bigger in cluster 2
      test_1$X3[i] <- rnorm(1,mean=2.0+0.8*test_1$X2[i],sd=sqrt(1))

      # test_1$R3_s[i] <- 2.5 - 0*test_1$X3[i] - 0.5*test_1$X2[i] +rnorm(1,mean=0,sd=1)
      # test_1$R[i] <- ifelse(test_1$R3_s[i]>0,1,0) #R=1 means observed.
      test_1$R[i] <- rbinom(1,1, pnorm(2.5 - 0*test_1$X3[i] - 0.5*test_1$X2[i]))
    }

    if (test_1$R[i] == 1) {test_1$X3_r[i] <- test_1$X3[i]}
  }
  # write.table(test_1$X3, file=paste("/Users/riceball/Desktop/PhD/Paper/Paper 1 - Mixture MNAR/temp/", "sample", t, "is being processed"))
  # table(test_1$cluster)
  # plot(density(test_1$X1),type="l", col="green", xlim=c(-3,6), ylim=c(0,0.8), main='')
  # par(new=TRUE)
  # plot(density(subset(test_1, test_1$cluster==1)$X1),type="l", col="red", xlim=c(-3,6), ylim=c(0,0.8), main='')
  # par(new=TRUE)
  # plot(density(subset(test_1, test_1$cluster==0)$X1),type="l", col="blue", xlim=c(-3,6), ylim=c(0,0.8), main='')
  # 
  # plot(density(test_1$X2),type="l", col="green", xlim=c(-3,6), ylim=c(0,0.8), main='')
  # par(new=TRUE)
  # plot(density(subset(test_1, test_1$cluster==1)$X2),type="l", col="red", xlim=c(-3,6), ylim=c(0,0.8), main='')
  # par(new=TRUE)
  # plot(density(subset(test_1, test_1$cluster==0)$X2),type="l", col="blue", xlim=c(-3,6), ylim=c(0,0.8), main='')
  # 
  # plot(density(test_1$X3),type="l", col="green", xlim=c(-3,6), ylim=c(0,0.8), main='')
  # par(new=TRUE)
  # plot(density(subset(test_1, test_1$cluster==1)$X3),type="l", col="red", xlim=c(-3,6), ylim=c(0,0.8), main='')
  # par(new=TRUE)
  # plot(density(subset(test_1, test_1$cluster==0)$X3),type="l", col="blue", xlim=c(-3,6), ylim=c(0,0.8), main='')
  # 
  # table(test_1$R)
  # summary(subset(test_1, test_1$cluster==1)$X1)
  # sd(subset(test_1, test_1$cluster==1)$X1)
  # summary(subset(test_1, test_1$cluster==0)$X1)
  # sd(subset(test_1, test_1$cluster==0)$X1)
  # lm(formula=X2~X1, subset(test_1, test_1$cluster==1))
  # lm(formula=X2~X1, subset(test_1, test_1$cluster==0))
  # lm(formula=X3~X2, subset(test_1, test_1$cluster==1))
  # lm(formula=X3~X2, subset(test_1, test_1$cluster==0))
  # glm(formula=R~X3+X2, subset(test_1, test_1$cluster==1), family = binomial(link = "probit"))
  # glm(formula=R~X3+X2, subset(test_1, test_1$cluster==0), family = binomial(link = "probit"))
  
  
  ######################################################################################################
  ######################################################################################################
  para <- data.frame(iteration=1:300,
                     prop=NA, #for binary cluster var "pik" (only consider two clusters)
                     #all below parameters are paired (for two clusters)
                     u00=NA, sd00=NA, u01=NA, sd01=NA, #for X1
                     b00=NA, b10=NA, sd0=NA, b01=NA, b11=NA, sd1=NA, #for X2|X1 (sd here is residual, not marginal sd)
                     beta00=NA, beta10=NA, sigma0=NA,  beta01=NA, beta11=NA, sigma1=NA, #for X3|X1,X2
                     a00=NA, a10=NA, a20=NA, a01=NA, a11=NA, a21=NA) #for missing model  (note, assume we already know the optimal form)
  #initial values
  m11 <- lm(formula=X2~X1, subset(test_1, test_1$cluster==1))$coefficient
  m12 <- lm(formula=X2~X1, subset(test_1, test_1$cluster==0))$coefficient
  m21 <- lm(formula=X3~X2, subset(test_1, test_1$cluster==1))$coefficient
  m22 <- lm(formula=X3~X2, subset(test_1, test_1$cluster==0))$coefficient
  m31 <- glm(formula=R~X3+X2, subset(test_1, test_1$cluster==1), family = binomial(link = "probit"))$coefficient
  # m32 <- glm(formula=R~X3+X2, subset(test_1, test_1$cluster==0), family = binomial(link = "probit"))$coefficient
  m32 <- glm(formula=R~X2, subset(test_1, test_1$cluster==0), family = binomial(link = "probit"))$coefficient #for MAR

  # para[1,2:24] <- c(table(test_1$cluster)[2]/n,
  #                   mean(subset(test_1, test_1$cluster==1)$X1),sd(subset(test_1, test_1$cluster==1)$X1), mean(subset(test_1, test_1$cluster==0)$X1),sd(subset(test_1, test_1$cluster==0)$X1),
  #                   m11[1], m11[2], 0.5, m12[1], m12[2], 0.5,
  #                   m21[1], m21[2], 1,   m22[1], m22[2], 1,
  #                   m31[1], m31[2],      m31[3], m32[1], m32[2], m32[3])
  para[1,2:24] <- c(table(test_1$cluster)[2]/n,
                    mean(subset(test_1, test_1$cluster==1)$X1),sd(subset(test_1, test_1$cluster==1)$X1), mean(subset(test_1, test_1$cluster==0)$X1),sd(subset(test_1, test_1$cluster==0)$X1),
                    m11[1], m11[2], 0.5,     m12[1], m12[2], 0.5,
                    m21[1], m21[2], 1,       m22[1], m22[2], 1,
                    m31[1], m31[2], m31[3],  m32[1], 0, m32[2])  #for MAR
  # 
  #note, although parameters for marginal Xobs part is directly identifiable, it is NOT directly identifiable under mixture
  
  #EM starts here:
  #E-step
  test_t <- test_1
  for (k in 1:200) { #iteration for EM
    
    if (k %in% c(5, 10, 20, 30, 40, 50, 80) | k %% 50 ==0) {print(paste("EM-iter", k, "is running"))}
    
    test_1 <- test_t[,1:7]
    
    #E-step 
    #1. compute E(Z=z|Xobs,R,current theta) 2. compute E(logf(Xmis|Xobs) & E(logP(R|X) w.r.p. to Xmis|Xobs,R)
    #to obtain compute above two Expectations, especially when the E has no closed-form in general, 
    #we use (gibbs) sampling to obtain latent X3 & Z in order to compute Q function, per subject and per current EM iteration's parameter value
    # start.time <- Sys.time()
    res <- foreach(i=1:n) %dopar% {
    # for (i in 1:n) {
      
      #for those with missing X3, we need to sample both X3 & Z
      if (test_1$R[i]==0){ #R=0 means missing
        
        #gibbs sampler for X3 & Z (cluster) or f(X3,Z|X1,X2,R;current parameter) per subject
        gibbs <- data.frame(iter=1:2000, X3=NA, Z=NA)
        gibbs[1,2:3] <- c(para$prop[k]*(para[k,13]+para[k,14]*test_1$X2[i])+(1-para$prop[k])*(para[k,16]+para[k,17]*test_1$X2[i]), 1)    #initial values
        
        for (index in 1:2000) {
          #1) sample Z (we calculate a binary proportion for z=1, since we only consider two clusters), 
          #per P(Z|Xmis_r-1, Xobs, R; current para) = P(Z)*P(X1|Z)*P(X2|X1,Z)*P(X3|X1,X2,Z)*P(R|X) / sum of ...
          #note, below prop1 is different from above prop, prop is a marginal para for mixing weight, while prop1 is per subject given X_i, R_i
          p1 <- para$prop[k] * #P(Z=1)
            dnorm(test_1$X1[i], mean=para$u00[k], sd=para$sd00[k]) * #f(X1|Z=1)
            dnorm(test_1$X2[i], mean=para[k,7]+para[k,8]*test_1$X1[i], sd=para$sd0[k]) *  #f(X2|X1,Z=1)
            dnorm(gibbs[index,2], mean=para[k,13]+para[k,14]*test_1$X2[i], sd=para$sigma0[k]) * #f(X3|X1,X2,Z=1)
            (1-pnorm(para[k,19]+para[k,20]*gibbs[index,2]+para[k,21]*test_1$X2[i]))  #P(R=0|X1,X2,X3,Z=1)
          # p2 <- (1-para$prop[k]) * 
          #   dnorm(test_1$X1[i], mean=para$u01[k], sd=para$sd01[k]) * 
          #   dnorm(test_1$X2[i], mean=para[k,10]+para[k,11]*test_1$X1[i], sd=para$sd1[k]) *
          #   dnorm(gibbs[index,2], mean=para[k,16]+para[k,17]*test_1$X2[i], sd=para$sigma1[k]) *
          #   (1-pnorm(para[k,22]+para[k,23]*gibbs[index,2]+para[k,24]*test_1$X2[i]))
          p2 <- (1-para$prop[k]) *   
            dnorm(test_1$X1[i], mean=para$u01[k], sd=para$sd01[k]) * 
            dnorm(test_1$X2[i], mean=para[k,10]+para[k,11]*test_1$X1[i], sd=para$sd1[k]) *
            dnorm(gibbs[index,2], mean=para[k,16]+para[k,17]*test_1$X2[i], sd=para$sigma1[k]) *
            (1-pnorm(para[k,22]+0*gibbs[index,2]+para[k,24]*test_1$X2[i]))    #for MAR
          prop1 <- p1/(p1+p2) #should be ~1 (1: labeled into 1st (true) cluster; 0: 2nd (wrong) cluster)
          
          #newly sampled Z_r for gibbs, then sample X3_r
          gibbs[index+1, 3] <- rbinom(1,1,prop1)
          
          #2) sample X3, we use ADS (since P(X3|X1,X2,R,Z_t) prop to P(X3|X1,X2,Z_r)*P(R|X1,X2,X3,Z_r))
          #note, the ADS here is conditioned on previously sampled Z_r (1 or 0), hence two separate ADS.
          if (gibbs[index+1, 3] == 1){
            #current EM iteration parameter value for cluster 1
            beta0 <- para[k,13]
            beta1 <- para[k,14]
            a0 <- para[k,19]
            a1 <- para[k,20]
            a2 <- para[k,21]
            sigma <- para$sigma0[k]
            
            #draw "1" sample for missing X3
            x2 <- test_1$X2[i]
            r <- test_1$R[i]
            f <- function(x3, beta0, beta1, sigma, tau, x2, a0, a1, a2, r){-(x3-beta0-beta1*x2)^2/(2*sigma^2) + 
                r*log(pnorm(a0+a1*x3+a2*x2))+(1-r)*log(1-pnorm(a0+a1*x3+a2*x2))} #log joint f 
            fprima <- function(x3, beta0, beta1, sigma, tau, x2, a0, a1, a2, r){-(x3-beta0-beta1*x2)/(sigma^2) + 
                sqrt(2/pi)*a1*r*exp(-(a0+a1*x3+a2*x2)^2/2)/(erf((a0+a1*x3+a2*x2)/sqrt(2))+1) - sqrt(2/pi)*a1*(r-1)*exp(-(a0+a1*x3+a2*x2)^2/2)/(erf((a0+a1*x3+a2*x2)/sqrt(2))-1) } #or not use erf but pnorm, same result.
            
            #two classic erros: 1. NA/inf error: log(1-pnorm) goes to -inf   2. trap, maximum updates: starting points are too few, need specify more
            cutoff <- (qnorm(1-10^-15)-a0-a2*x2)/a1
            cutoff2 <- (qnorm(10^-15)-a0-a2*x2)/a1
            lb <- -10   #original bound, could choose a relatively large number (to avoid trap error), since it will be further protected by cutoff to avoid NA/inf error.
            ub <- 10  #original bound
            if (a1<0) {
              range <- c(seq(max(ceiling(cutoff),lb),(max(ceiling(cutoff),lb)+min(floor(cutoff2),ub))/2, by=0.5), 
                         rev(seq(min(floor(cutoff2),ub), (max(ceiling(cutoff),lb)+min(floor(cutoff2),ub))/2, by=-0.5)))
              # range <- c(seq(max(cutoff,lb), 0, by=0.5), rev(seq(min(cutoff2,ub),0, by=-0.5)))
              mysample <- ars(n=1,f,fprima, m=length(range), x=range, lb=TRUE,xlb=cutoff,ub=TRUE,xub=cutoff2,
                              beta0=beta0, beta1=beta1, sigma=sigma, x2=x2, a0=a0, a1=a1, a2=a2, r=r)
            }else if (a1>0) {
              range <- c(seq(max(ceiling(cutoff2),lb),(max(ceiling(cutoff2),lb)+min(floor(cutoff),ub))/2, by=0.5), 
                         rev(seq(min(floor(cutoff),ub), (max(ceiling(cutoff2),lb)+min(floor(cutoff),ub))/2, by=-0.5)))
              # range <- c(seq(max(cutoff2,lb), 0, by=0.5), rev(seq(min(cutoff,ub),0, by=-0.5)))
              mysample <- ars(n=1,f,fprima, m=length(range), x=range, lb=TRUE,xlb=cutoff2,ub=TRUE,xub=cutoff,
                              beta0=beta0, beta1=beta1, sigma=sigma, x2=x2, a0=a0, a1=a1, a2=a2, r=r)
            }
            #newly sampled X3_r
            gibbs[index+1, 2] <- mysample
            remove(mysample)
            
          }else if (gibbs[index+1, 3] == 0) {
            #current EM iteration parameter value for cluster 2
            # beta0 <- para[k,16]
            # beta1 <- para[k,17]
            # a0 <- para[k,22]
            # a1 <- 0.00001  #note, since we hypothesize cluster 2 to a MAR, another simplied method is that we just a1 ==0 or say 0.0001, 
            # #but previous simu shows using a more saturated MNAR model has as good result as directly using a MAR model.
            # a2 <- para[k,24]
            # sigma <- para$sigma1[k]
            
            #draw "1" sample for missing X3
            x2 <- test_1$X2[i]
            # r <- test_1$R[i]
            # f <- function(x3, beta0, beta1, sigma, tau, x2, a0, a1, a2, r){-(x3-beta0-beta1*x2)^2/(2*sigma^2) + 
            #     r*log(pnorm(a0+a1*x3+a2*x2))+(1-r)*log(1-pnorm(a0+a1*x3+a2*x2))} #log joint f 
            # fprima <- function(x3, beta0, beta1, sigma, tau, x2, a0, a1, a2, r){-(x3-beta0-beta1*x2)/(sigma^2) + 
            #     sqrt(2/pi)*a1*r*exp(-(a0+a1*x3+a2*x2)^2/2)/(erf((a0+a1*x3+a2*x2)/sqrt(2))+1) - sqrt(2/pi)*a1*(r-1)*exp(-(a0+a1*x3+a2*x2)^2/2)/(erf((a0+a1*x3+a2*x2)/sqrt(2))-1) } #or not use erf but pnorm, same result.
            # 
            # #two classic erros: 1. NA/inf error: log(1-pnorm) goes to -inf   2. trap, maximum updates: starting points are too few, need specify more
            # cutoff <- (qnorm(1-10^-15)-a0-a2*x2)/a1
            # cutoff2 <- (qnorm(10^-15)-a0-a2*x2)/a1
            # lb <- -10   #original bound, could choose a relatively large number (to avoid trap error), since it will be further protected by cutoff to avoid NA/inf error.
            # ub <- 10  #original bound
            # if (a1<0) {
            #   range <- c(seq(max(ceiling(cutoff),lb),(max(ceiling(cutoff),lb)+min(floor(cutoff2),ub))/2, by=0.5), 
            #              rev(seq(min(floor(cutoff2),ub), (max(ceiling(cutoff),lb)+min(floor(cutoff2),ub))/2, by=-0.5)))
            #   # range <- c(seq(max(cutoff,lb), 0, by=0.5), rev(seq(min(cutoff2,ub),0, by=-0.5)))
            #   mysample <- ars(n=1,f,fprima, m=length(range), x=range, lb=TRUE,xlb=cutoff,ub=TRUE,xub=cutoff2,
            #                   beta0=beta0, beta1=beta1, sigma=sigma, x2=x2, a0=a0, a1=a1, a2=a2, r=r)
            # }else if (a1>0) {
            #   range <- c(seq(max(ceiling(cutoff2),lb),(max(ceiling(cutoff2),lb)+min(floor(cutoff),ub))/2, by=0.5), 
            #              rev(seq(min(floor(cutoff),ub), (max(ceiling(cutoff2),lb)+min(floor(cutoff),ub))/2, by=-0.5)))
            #   # range <- c(seq(max(cutoff2,lb), 0, by=0.5), rev(seq(min(cutoff,ub),0, by=-0.5)))
            #   mysample <- ars(n=1,f,fprima, m=length(range), x=range, lb=TRUE,xlb=cutoff2,ub=TRUE,xub=cutoff,
            #                   beta0=beta0, beta1=beta1, sigma=sigma, x2=x2, a0=a0, a1=a1, a2=a2, r=r)
            # }
            
            #for MAR: note, in this case, P(X3|X1,X2,R,Z) just prop to P(X3|X1,X2,Z), no ARS needed.
            #newly sampled X3_r
            gibbs[index+1, 2] <- rnorm(1, mean=para[k,16]+para[k,17]*test_1$X2[i], sd=para$sigma1[k])
          }
        }
        
        #obtain (Z,X3)|X1,X2,R, current para, for E-step
        #note, to compute the parenthesis part expectation (w.r.p to X3|Z,Xobs), we will use two different subsamples (per Z) from above marginal sample of X3
        #i.e., for M-step, each subject will be used twice in two clusters, each with a double-weight (=marginal pi_k*weight from subsample) 
        ss <- seq(1001,2000,by=2)  #burn-in=1000, thinning=2
        
        test_1[i,8] <- sum(gibbs$Z[ss] == 1)/500  #for P(Z|Xobs,R;current para), 1st weight
        test_1[i,9] <- sum(gibbs$Z[ss] == 1)
        colnames(test_1)[8] <- c("prop1")
        colnames(test_1)[9] <- c("n1")
        
        sub1 <- subset(gibbs$X3, gibbs$Z==1 & gibbs$iter %in% ss) #for P(Xmis|Xobs,R,Z;current para), its sample size =1/2nd weight
        sub2 <- subset(gibbs$X3, gibbs$Z==0 & gibbs$iter %in% ss) #for P(Xmis|Xobs,R,Z;current para), its sample size =1/2nd weight
        if (length(sub1)>0 & length(sub2)>0){
          test_1[i,10:(8+length(sub1)+1)] <- sub1
          test_1[i,(8+length(sub1)+2):(8+length(sub1)+1+length(sub2))] <- sub2
        }else if (length(sub1)==0 & length(sub2)>0){
          test_1[i,(8+length(sub1)+2):(8+length(sub1)+1+length(sub2))] <- sub2
        }else if (length(sub2)==0 & length(sub1)>0){
          test_1[i,10:(8+length(sub1)+1)] <- sub1
        }
        remove(sub1)
        remove(sub2)
        remove(gibbs)
        #note, in this regard, every subj has different N for two subsamples
        
        
      }else if (test_1$R[i]==1){ #for those without missing X3, we only need to consider Z (and no gibbs)
        
        #calculate exact P(Z|X,R;current para) (only difference with above gibbs is that now Xmis is fixed, no gibbs on Xmis)
        #P(Z|Xmis, Xobs, R; current para) = P(Z)*P(X1|Z)*P(X2|X1,Z)*P(X3|X1,X2,Z)*P(R|X) / sum of ...
        #note, for this P(Z|X,R;current para), when marginal P(Z=1) is ~1, P(Z=1|X,R;current para) will almost surely dominate by term P(Z=1), 
        #i.e., almost all data (even for Z=0) will be labeled into cluster 1 just because of this term P(Z=1) >> P(Z=0), 
        #this phenomnon will be attenuated only when two clusters have extremely different P(X|Z) and/or P(R|X,Z) to reverse this term's effect.
        #compute P(X3,X1,X2,Z=z,R;current para) for Z=1,0; then divided by sum or P(X3,X1,X2,R;current para).
        p1 <- para$prop[k] * #P(Z=1)
          dnorm(test_1$X1[i], mean=para$u00[k], sd=para$sd00[k]) * #f(X1|Z=1)
          dnorm(test_1$X2[i], mean=para[k,7]+para[k,8]*test_1$X1[i], sd=para$sd0[k]) *  #f(X2|X1,Z=1)
          dnorm(test_1$X3_r[i], mean=para[k,13]+para[k,14]*test_1$X2[i], sd=para$sigma0[k]) * #f(X3|X1,X2,Z=1)
          pnorm(para[k,19]+para[k,20]*test_1$X3_r[i]+para[k,21]*test_1$X2[i])  #P(R=1|X1,X2,X3,Z=1)
        # p2 <- (1-para$prop[k]) * 
        #   dnorm(test_1$X1[i], mean=para$u01[k], sd=para$sd01[k]) * 
        #   dnorm(test_1$X2[i], mean=para[k,10]+para[k,11]*test_1$X1[i], sd=para$sd1[k]) *
        #   dnorm(test_1$X3_r[i], mean=para[k,16]+para[k,17]*test_1$X2[i], sd=para$sigma1[k]) *
        #   pnorm(para[k,22]+para[k,23]*test_1$X3_r[i]+para[k,24]*test_1$X2[i])
        p2 <- (1-para$prop[k]) * 
          dnorm(test_1$X1[i], mean=para$u01[k], sd=para$sd01[k]) * 
          dnorm(test_1$X2[i], mean=para[k,10]+para[k,11]*test_1$X1[i], sd=para$sd1[k]) *
          dnorm(test_1$X3_r[i], mean=para[k,16]+para[k,17]*test_1$X2[i], sd=para$sigma1[k]) *
          pnorm(para[k,22]+0*test_1$X3_r[i]+para[k,24]*test_1$X2[i])   #for MAR
        
        #obtain P(Z=1|Xmis,Xobs,R;current para), 1st (and sole) weight 
        test_1[i,8] <-  p1/(p1+p2) 
        colnames(test_1)[8] <- c("prop1")
        
      }
      list(test_1[i,])
      
    } #complete E-step 
  
    # test_1 <- ldply (res, data.frame) #same effect as below syntax, but slower
    list_data <- Map(as.data.frame, res) 
    test_1  <- rbindlist(list_data, fill=TRUE)
    remove(list_data)
    remove(res)
    # end.time <- Sys.time()
    # time.taken <- end.time - start.time
    # time.taken
    
    #compare sampled X3|X2,para vs. true X3|X2,para 
    #注意，其实这个比较是不准确的，因为gibbs（ARS）得到的是X3|X2,R，而下面的true只是X3|X2.
    #这两个只有当所有X3 + observed X2对应的P(R|.)都几乎一样时，才会近似相等。
    #(比如对于i=7，summary(1-pnorm(a0+a1*sub1+a2*x2))其实不都一样，所以两个不一样；而对于i=22，summary(1-pnorm(para[k,16]+para[k,17]*sub1))几乎都是很小的值，所以两者几乎一样）
    # t <- 24 #subj ID
    # if (test_1$cluster[t] ==0){
    #   plot(density(t(test_1[t,(8+test_1$n1[t]+2):(8+test_1$n1[t]+1+500-test_1$n1[t])])),col='red',xlim=c(-2,6), ylim=c(0,1), main='') #sampled d
    #   par(new=TRUE)
    #   lines(density(rnorm(10000,mean=para[k,16]+para[k,17]*test_1$X2[t], sd=para$sigma1[k])), col='blue') #per current para
    #   par(new=TRUE)
    #   lines(density(rnorm(10000,mean=para[1,16]+para[1,17]*test_1$X2[t], sd=para$sigma1[1])), col='green') #true d
    # }else if (test_1$cluster[t] ==1){
    #   plot(density(t(test_1[t,10:(8+test_1$n1[t]+1)])),col='red',xlim=c(-2,6), ylim=c(0,1), main='') #sampled d
    #   par(new=TRUE)
    #   lines(density(rnorm(10000,mean=para[k,13]+para[k,14]*test_1$X2[t], sd=para$sigma0[k])), col='blue') #per current para
    #   par(new=TRUE)
    #   lines(density(rnorm(10000,mean=para[1,13]+para[1,14]*test_1$X2[t], sd=para$sigma0[1])), col='green') #true d
    #   #summary(t(test_1[t,12:(11+test_1$n1[t])]))
    # }

    #M-step
    #parameters to be updated: 
    #1) pi_K  
    #2) u & sd for X1  
    #3) b0 b1 & sd for X2|X1   
    #4) beta0 beta1 & sigma for X3|X2
    #5) a0 a1 a2 for R|X2,X3
    # para <- data.frame(iteration=1:1500, 
    #                    prop=NA, #for binary cluster var "pik" (only consider two clusters)
    #                    #all below parameters are paired (for two clusters)
    #                    u00=NA, sd00=NA, u01=NA, sd01=NA, #for X1 
    #                    b00=NA, b10=NA, sd0=NA, b01=NA, b11=NA, sd1=NA, #for X2|X1 (sd here is residual, not marginal sd)
    #                    beta00=NA, beta10=NA, sigma0=NA,  beta01=NA, beta11=NA, sigma1=NA, #for X3|X1,X2  
    #                    a00=NA, a10=NA, a20=NA, a01=NA, a11=NA, a21=NA) #for missing model  (note, assume we already know the optimal form)
    
    #re-format the data into long-form
    long11 <- long12 <- subset(test_1, test_1$R==1, select=c(ID, X3_r,X2, X1, R, prop1)) #CC
    long11$cluster <- 1
    long12$prop1 <- 1- long12$prop1
    long12$cluster <- 0
    long1 <- rbind(long11, long12) 
    long1 <- long1[order(long1$ID),]
    long1$w <- 1 
    long1 <- long1[,-1]
    remove(long11)
    remove(long12)
    
    long2 <- data.frame(ID=1:(nrow(subset(test_1,test_1$R==0))*500),X3_r=NA, X2=NA, X1=NA, R=0, cluster=NA, prop1=NA,w=NA) #with missing X3
    temp <- subset(test_1, test_1$R==0, select=-c(ID,X3,X3_r,X1,X2,R, prop1, cluster,n1))
    long2$X3_r <- as.vector(t(temp))
    for (i in 1:nrow(subset(test_1,test_1$R==0))){
      long2$X1[((i-1)*500+1):((i-1)*500+1+500-1)] <- subset(test_1, test_1$R==0, select=X1)$X1[i]
      long2$X2[((i-1)*500+1):((i-1)*500+1+500-1)] <- subset(test_1, test_1$R==0, select=X2)$X2[i]
      if (subset(test_1,test_1$R==0)$n1[i] ==0){  #all sampled Z are in cluster 0
        long2$prop1[(((i-1)*500+1)):((i-1)*500+500)]  <- 1
        long2$w[(((i-1)*500+1)):((i-1)*500+500)]  <- 1/500
        long2$cluster[(((i-1)*500+1)):((i-1)*500+500)]  <- 0
      }else if (subset(test_1,test_1$R==0)$n1[i] ==500){ #all sampled Z are in cluster 1
        long2$prop1[((i-1)*500+1):(((i-1)*500+1)+500-1)]  <-  1
        long2$w[((i-1)*500+1):(((i-1)*500+1)+500-1)]  <-  1/500
        long2$cluster[((i-1)*500+1):(((i-1)*500+1)+500-1)]  <-  1
      }else {
        long2$prop1[((i-1)*500+1):(((i-1)*500+1)+subset(test_1,test_1$R==0)$n1[i]-1)]  <-  subset(test_1,test_1$R==0)$prop1[i]
        long2$prop1[(((i-1)*500+1)+subset(test_1,test_1$R==0)$n1[i]):((i-1)*500+500)]  <- 1-subset(test_1,test_1$R==0)$prop1[i]
        long2$w[((i-1)*500+1):(((i-1)*500+1)+subset(test_1,test_1$R==0)$n1[i]-1)]  <-  1/subset(test_1,test_1$R==0)$n1[i]
        long2$w[(((i-1)*500+1)+subset(test_1,test_1$R==0)$n1[i]):((i-1)*500+500)]  <- 1/(500-subset(test_1,test_1$R==0)$n1[i])
        long2$cluster[((i-1)*500+1):(((i-1)*500+1)+subset(test_1,test_1$R==0)$n1[i]-1)]  <-  1
        long2$cluster[(((i-1)*500+1)+subset(test_1,test_1$R==0)$n1[i]):((i-1)*500+500)]  <-  0
      }
    }
    long2 <- long2[,-1]
    remove(temp)
    
    long <- rbind(long1, long2) 
    long$w2 <- long$w*long$prop1
    for (i in 1:length(long$X3_r)){
      if (is.na(long$w2[i])=='TRUE') {long$w2[i] <- 0}
    }
    # sum(long$w2)
    
    #1) update marginal pi_k
    #note, for pi_k, only the 1st weight is needed (E(I(Zi=k)) w.r.p to Z|Xobs,R), and it has a closed-form MLE solution.
    #for those without missing, the weight is exact; for those with missing, the weight is computed via gibbs sampler.
    para$prop[k+1] <- sum(test_1$prop1)/n   #n=(sum(test_1$prop1)+sum(1-test_1$prop1))
    
    #2) update u & sd for X1  
    #note, here we will use two weights (one from E(I(Zi=k)), one from gibbs sampler for X3|X2,current para)
    #2.1 update X1's para (has closed-form)
    #u00=NA, sd00=NA, u01=NA, sd01=NA, for X1 
    d1 <- subset(long,long$cluster==1)
    d2 <- subset(long,long$cluster==0)
    # para$u00[k+1] <- sum(d1$X1*d1$w2)/sum(d1$w2)
    # para$u01[k+1] <- sum(d2$X1*d2$w2)/sum(d2$w2)
    fit_X11 <- lm(formula=X1~1, data=d1, weights=w2)
    fit_X12 <- lm(formula=X1~1, data=d2, weights=w2)
    d1$fitted1 <- fit_X11$fitted.values
    d2$fitted1 <- fit_X12$fitted.values
    
    para$u00[k+1] <- fit_X11[["coefficients"]]
    para$u01[k+1] <- fit_X12[["coefficients"]]
    
    para$sd00[k+1] <- sqrt(sum((d1$X1-d1$fitted1)^2*d1$w2) /sum(d1$w2))
    para$sd01[k+1] <- sqrt(sum((d2$X1-d2$fitted1)^2*d2$w2) /sum(d2$w2))
    
    #3) update b & sd for X2|X1
    #b00=NA, b10=NA, sd0=NA, b01=NA, b11=NA, sd1=NA, #for X2|X1 (sd here is residual, not marginal sd)
    fit_X21 <- lm(formula=X2~X1, data=d1, weights=w2)
    fit_X22 <- lm(formula=X2~X1, data=d2, weights=w2)
    d1$fitted2 <- fit_X21$fitted.values
    d2$fitted2 <- fit_X22$fitted.values
    
    para$b00[k+1] <- fit_X21[["coefficients"]][1]
    para$b10[k+1] <- fit_X21[["coefficients"]][2]
    para$b01[k+1] <- fit_X22[["coefficients"]][1]
    para$b11[k+1] <- fit_X22[["coefficients"]][2]
    
    para$sd0[k+1] <- sqrt(sum((d1$X2-d1$fitted2)^2*d1$w2) /sum(d1$w2))
    para$sd1[k+1] <- sqrt(sum((d2$X2-d2$fitted2)^2*d2$w2) /sum(d2$w2))
    
    #4) beta0 beta1 & sigma for X3|X2
    #beta00=NA, beta10=NA, sigma0=NA,  beta01=NA, beta11=NA, sigma1=NA, #for X3|X1,X2  
    # lm(formula=X3~X2, data=subset(test_1,test_1$cluster==1)) true
    # lm(formula=X3~X2, data=subset(test_1,test_1$cluster==0)) true
    fit_X31 <- lm(formula=X3_r~X2, data=d1, weights=w2)
    fit_X32 <- lm(formula=X3_r~X2, data=d2, weights=w2)
    d1$fitted3 <- fit_X31$fitted.values
    d2$fitted3 <- fit_X32$fitted.values
    
    para$beta00[k+1] <- fit_X31[["coefficients"]][1]
    para$beta10[k+1] <- fit_X31[["coefficients"]][2]
    para$beta01[k+1] <- fit_X32[["coefficients"]][1]
    para$beta11[k+1] <- fit_X32[["coefficients"]][2]
    
    para$sigma0[k+1] <- sqrt(sum((d1$X3_r-d1$fitted3)^2*d1$w2) /sum(d1$w2))
    para$sigma1[k+1] <- sqrt(sum((d2$X3_r-d2$fitted3)^2*d2$w2) /sum(d2$w2))
    
    #5) a0 a1 a2 for R|X2,X3
    #a00=NA, a10=NA, a20=NA, a01=NA, a11=NA, a21=NA) #for missing model  (note, assume we already know the optimal form)
    # glm(formula=R~X3+X2, subset(test_1,test_1$cluster==1), family = binomial(link = "probit"))
    # glm(formula=R~X3+X2, subset(test_1,test_1$cluster==0), family = binomial(link = "probit"))
    fit_R1 <- glm(formula=R~X3_r+X2, d1, family = binomial(link = "probit"), weights=w2)
    # fit_R2 <- glm(formula=R~X3_r+X2, d2, family = binomial(link = "probit"), weights=w2)
    fit_R2 <- glm(formula=R~X2, d2, family = binomial(link = "probit"), weights=w2)  #for MAR
    
    para$a00[k+1] <- fit_R1[["coefficients"]][1]
    para$a10[k+1] <- fit_R1[["coefficients"]][2]
    para$a20[k+1] <- fit_R1[["coefficients"]][3]
    para$a01[k+1] <- fit_R2[["coefficients"]][1]
    para$a11[k+1] <- 0
    para$a21[k+1] <- fit_R2[["coefficients"]][2]
    
    # #likelihood:
    # sum(d2$w2*(d2$R*log(pnorm(para$a01[k]+para$a11[k]*d2$X3_r+para$a21[k]*d2$X2)) +
    #            (1-d2$R)*log(1-pnorm(para$a01[k]+para$a11[k]*d2$X3_r+para$a21[k]*d2$X2))))
    # sum(d2$w2*(d2$R*log(pnorm(para$a01[1]+para$a11[1]*d2$X3_r+para$a21[1]*d2$X2)) +
    #              (1-d2$R)*log(1-pnorm(para$a01[1]+para$a11[1]*d2$X3_r+para$a21[1]*d2$X2))))
    # 
    # sum(d1$w2*(d1$R*log(pnorm(para$a00[k]+para$a10[k]*d1$X3_r+para$a20[k]*d1$X2)) +
    #              (1-d1$R)*log(1-pnorm(para$a00[k]+para$a10[k]*d1$X3_r+para$a20[k]*d1$X2))))
    # sum(d1$w2*(d1$R*log(pnorm(para$a00[1]+para$a10[1]*d1$X3_r+para$a20[1]*d1$X2)) +
    #              (1-d1$R)*log(1-pnorm(para$a00[1]+para$a10[1]*d1$X3_r+para$a20[1]*d1$X2))))
    # 
    # sum(1*(s2$R*log(pnorm(para$a01[k]+para$a11[k]*s2$X3+para$a21[k]*s2$X2)) +
    #              (1-s2$R)*log(1-pnorm(para$a01[k]+para$a11[k]*s2$X3+para$a21[k]*s2$X2))))
    # sum(1*(s2$R*log(pnorm(para$a01[1]+para$a11[1]*s2$X3+para$a21[1]*s2$X2)) +
    #              (1-s2$R)*log(1-pnorm(para$a01[1]+para$a11[1]*s2$X3+para$a21[1]*s2$X2))))
    
    remove(d1)
    remove(d2)
    remove(long)
    remove(long1)
    remove(long2)
    
    remove(fit_X11)
    remove(fit_X12)
    remove(fit_X21)
    remove(fit_X22)
    remove(fit_X31)
    remove(fit_X32)
    remove(fit_R1)
    remove(fit_R2)
    
    #convergence
    #squared distance for all parameters are <10^-5 from previous iteration (equivalent to say the key 17 para has a 0.0002 change)
    if (k>=40) {
      if (sum((para[k,2:18]-para[k-1,2:18])^2) <10^(-6)|
          sum((para[k,-1]-para[k-1,-1])^2) <10^(-6) |
          (k==200))  break
    }
    # para <- paraa[[9]]
    
  } #one EM iteration completes. 
  
  testt[[t]] <- test_1
  paraa[[t]] <- subset(para, is.na(para$prop)==0)
  
  # list(test_1, subset(para, is.na(para$prop)==0))
  
}
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken


###############################################################################################################################
#average result
N <- t-1
N <- 500
summary <- data.frame(sample=1:500, 
                      prop=NA, #for binary cluster var "pik" (only consider two clusters)
                      #all below parameters are paired (for two clusters)
                      u00=NA, sd00=NA, u01=NA, sd01=NA, #for X1 
                      b00=NA, b10=NA, sd0=NA, b01=NA, b11=NA, sd1=NA, #for X2|X1 (sd here is residual, not marginal sd)
                      beta00=NA, beta10=NA, sigma0=NA,  beta01=NA, beta11=NA, sigma1=NA, #for X3|X1,X2  
                      a00=NA, a10=NA, a20=NA, a01=NA, a11=NA, a21=NA) #for missing model  (note, assume we already know the optimal form)
for (i in 1:N){
  s <- subset(paraa[[i]], is.na(paraa[[i]]$prop)==0)
  summary[i,2:24] <- s[nrow(s),2:24]
  summary[i,25] <- nrow(s)
}
summary(summary)


truesummary <- data.frame(sample=1:500, 
                          prop=NA, #for binary cluster var "pik" (only consider two clusters)
                          #all below parameters are paired (for two clusters)
                          u00=NA, sd00=NA, u01=NA, sd01=NA, #for X1 
                          b00=NA, b10=NA, sd0=NA, b01=NA, b11=NA, sd1=NA, #for X2|X1 (sd here is residual, not marginal sd)
                          beta00=NA, beta10=NA, sigma0=NA,  beta01=NA, beta11=NA, sigma1=NA, #for X3|X1,X2  
                          a00=NA, a10=NA, a20=NA, a01=NA, a11=NA, a21=NA) #for missing model  (note, assume we already know the optimal form)
for (i in 1:N){
  m11 <- lm(formula=X2~X1, subset(testt[[i]], testt[[i]]$cluster==1))$coefficient
  m12 <- lm(formula=X2~X1, subset(testt[[i]], testt[[i]]$cluster==0))$coefficient
  m21 <- lm(formula=X3~X2, subset(testt[[i]], testt[[i]]$cluster==1))$coefficient
  m22 <- lm(formula=X3~X2, subset(testt[[i]], testt[[i]]$cluster==0))$coefficient
  m31 <- glm(formula=R~X3+X2, subset(testt[[i]], testt[[i]]$cluster==1), family = binomial(link = "probit"))$coefficient
  # m32 <- glm(formula=R~X3+X2, subset(testt[[i]], testt[[i]]$cluster==0), family = binomial(link = "probit"))$coefficient
  m32 <- glm(formula=R~X2, subset(testt[[i]], testt[[i]]$cluster==0), family = binomial(link = "probit"))$coefficient #for MAR
  
  truesummary[i,2] <- table(testt[[i]]$cluster)[2]/n
  truesummary[i,3] <- mean(subset(testt[[i]], testt[[i]]$cluster==1)$X1)
  truesummary[i,4] <- sd(subset(testt[[i]], testt[[i]]$cluster==1)$X1)
  truesummary[i,5] <- mean(subset(testt[[i]], testt[[i]]$cluster==0)$X1)
  truesummary[i,6] <- sd(subset(testt[[i]], testt[[i]]$cluster==0)$X1)
  truesummary[i,7:8] <- m11
  truesummary[i,10:11] <- m12
  truesummary[i,13:14] <- m21
  truesummary[i,16:17] <- m22
  truesummary[i,19:21] <- m31
  truesummary[i,22] <- m32[1]
  truesummary[i,24] <- m32[2]
}
summary(truesummary)


s_1 <- c(2,3,4,7,8,9,13,14,15,19,20,21)
s_2 <- c(2,5,6,10,11,12,16,17,18,22,23,24)

#bias, SE, MSE
true1 <- c(prop, 0,0.5, 0.8,0.5,0.5, 0.8,0.8,1, 1,0.8,-0.5) #10 miss
# true1 <- c(prop, 0,0.5, 0.8,0.5,0.5, 0.8,0.8,1,   2.5,-0.9,-0.5) #30 miss
true2 <- c(1-prop, 1.5,0.5, 2,0.5,0.5, 2,0.8,1,   2.5, 0,-0.5)
temp <- summary[nrow(summary),s_1]
temp[1,] <- mapply(mean,summary[,s_1],na.rm=T)
temp[2,] <- temp[1,] - true1
temp[3,] <- apply(summary[1:N,s_1], 2, sd) #SE
temp[4,] <- temp[2,]^2+temp[3,]^2

temp[5,] <- mapply(mean,summary[,s_2],na.rm=T)
temp[5,1] <- 1-temp[5,1]
temp[6,] <- temp[5,] - true2
temp[7,] <- apply(summary[1:N,s_2], 2, sd) #SE
temp[8,] <- temp[6,]^2+temp[7,]^2

temp <- round(temp,4)

#boxplot
summaryy <- summary[,s_1]
for (i in 1:N){
  summaryy[i,] <- summaryy[i,] - true1
}
summaryy2 <- summary[,s_2]
summaryy2$prop <- 1-summaryy2$prop 
for (i in 1:N){
  summaryy2[i,] <- summaryy2[i,] - true2
}
options(scipen=10000)
par(mfrow = c(1, 2))
#outcome model part
boxplot(summaryy[,1],summaryy[,2],summaryy[,3],summaryy[,4],summaryy[,5],summaryy[,6],
        summaryy2[,1],summaryy2[,2],summaryy2[,3],summaryy2[,4],summaryy2[,5],summaryy2[,6],
        at =c(1, 3,5, 7,9,11, 
              2, 4,6, 8,10,12), outline=FALSE,
        col=c("white","white","white","white","white","white",
              "grey","grey","grey","grey","grey","grey"),
        ylab ="",
        ylim = c(-0.6, 0.6),
        # main="Parameters related to X1 and X2 (mixture setting 1)",
        names=c("prop",expression(mu[1]),expression(sigma[1]),expression(alpha[0]),expression(alpha[1]),expression(sigma[2]),
              "","","","","",""))
abline(h=0)
points(x=c(1, 3,5, 7,9,11,    
           2, 4,6, 8,10,12), 
       y=c(temp[2,1], temp[2,c(2,3)], temp[2,c(4,5,6)], 
           temp[6,1], temp[6,c(2,3)], temp[6,c(4,5,6)]), col = "red")
# text(x=c(1, 3,5, 7,9,11,    13,15,17,
#          2, 4,6, 8,10,12,   14,16,18), 
#      # cex = 0.8,
#      y=c(round(temp[2,1:9],4)+0.3,round(temp[6,1:9],4)+0.3), labels = c(round(temp[2,1:9],4),round(temp[6,1:9],4)))

#selection model part
boxplot(summaryy[,7],summaryy[,8],summaryy[,9], summaryy[,10],summaryy[,11],summaryy[,12],
        summaryy2[,7],summaryy2[,8],summaryy2[,9], summaryy2[,10],summaryy2[,11],summaryy2[,12],
        at =c(1, 3,5, 7,9,11,    
              2, 4,6, 8,10,12), outline=FALSE,
        col=c("white","white","white","white","white","white",
              "grey","grey","grey","grey","grey","grey"),
        ylab ="",
        ylim = c(-3.8, 3.8),
        # main="Parameters related to X3 (mixture setting 1)",
        names=c(expression(beta[0]),expression(beta[1]),expression(sigma[3]), expression(phi[0]),expression(phi[1]),expression(phi[2]),
                "","","","","",""))
abline(h=0)
points(x=c(1, 3,5, 7,9,11,    
           2, 4,6, 8,10,12), 
       y=c(temp[2,c(7,8,9)], temp[2,c(10,11,12)], 
           temp[6,c(7,8,9)], temp[6,c(10,11,12)]), col = "red")
# text(x=c(1, 3,5,
#          2, 4,6), 
#      # cex = 0.8,
#      y=c(round(temp[2,10:12],4)+0.3,round(temp[6,10:12],4)+0.3), labels = c(round(temp[2,10:12],4),round(temp[6,10:12],4)))


write.csv(temp, "/Users/riceball/Desktop/PhD/Paper/Paper 1 - Mixture MNAR/temp.csv", row.names=FALSE)

#for latex output
noquote(paste(format(round(as.numeric(cbind(temp[2,])),4)), "&"))
noquote(paste(format(round(as.numeric(cbind(temp[3,])),4)), "&"))
noquote(paste(format(round(as.numeric(cbind(temp[4,])),4)), "&"))

noquote(paste(format(round(as.numeric(cbind(temp[6,])),4)), "&"))
noquote(paste(format(round(as.numeric(cbind(temp[7,])),4)), "&"))
noquote(paste(format(round(as.numeric(cbind(temp[8,])),4)), "&"))


#clustering 
c1 <- rep(NA)
c2 <- rep(NA)
for (t in 1:N){
  c1[t] <- 1-length(subset(testt[[t]],testt[[t]]$cluster==1 & testt[[t]]$prop1<0.5)$ID)/sum(testt[[t]]$cluster)
  c2[t] <- 1-length(subset(testt[[t]],testt[[t]]$cluster==0 & testt[[t]]$prop1>0.5)$ID)/(n-sum(testt[[t]]$cluster))
}
summary(c1)
summary(c2)

####################################################################################################################################
#impute Xmis
for (p in 1:20){
  name <- paste("test_imp", p, sep = "_")
  tempp <- assign(name, test_1[,1:7])
  
  for (i in 1:500){
  
    
  }
 
  
  remove(tempp)
}

remove(test_imp_1)










