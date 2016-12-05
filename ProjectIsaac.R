# load file
strike <- read.csv("Documents/Waterloo/4A/331/project/strikes_clean.csv")

# Taken from Will

# encode categoricals as Factors
strike$Centr <- as.factor(strike$Centr)
strike$Country <- as.factor(strike$Country)
strike$Demo <- log(strike$Demo)
strike$Dens <- log(strike$Dens)

#strike$Unemp[strike$Unemp<=0] <- 0.00000000001
#strike$Infl[strike$Infl<=0] <- 0.00000000001
#strike$Demo[strike$Demo<=0] <- 0.00000000001
#strike$Dens[strike$Dens<=0] <- 0.00000000001

correlMatrix <- cbind(strike$Unemp, strike$Infl, strike$Demo, strike$Dens, strike$Centr)
cor(correlMatrix)

# Reason that these should affect strike
M1 = lm(Strike ~ Centr + Infl + Centr:Infl, data = strike)
# anova shows unemp doesn't add much 
anova(M1)
# remove Unemp from M1 to get current model that works well 
#M2 = lm(Strike ~ Centr + Infl + log(Dens) + log(Demo) + Centr:Infl, data = strike)
# M2 = lm(Strike ~ Centr + Infl + log(Dens) + log(Demo) + Centr:Infl + Centr:log(Dens) + Infl:log(Demo), data = strike)
# basic model
M0 <- lm(Strike ~ 1, data = strike)
# full model without country - interaction effects

Mfull <- lm(Strike ~ (. - Country - Year)^2, data = strike)
#Mfull <- lm(Strike ~ (. - Country - Year)^2, data = strike)

# forward step model from intercept to full
Mfwd <- step(M0, scope=list(lower=M0, upper=Mfull), direction='forward', trace=FALSE)
# backward step model from full to intercept
Mback <- step(Mfull, scope=list(lower=M0, upper=Mfull), direction='backward', trace=FALSE)
# step model starting from model grenerated from exploration
Mstep <- step(M1, scope=list(lower=M0, upper=Mfull), direction='both', trace=FALSE)


models <- list(M0, M1, M2, Mfwd, Mback, Mstep)

# compare M1 and M2

M1 <- M1
M2 <- M2
Mnames <- expression(M[BACK], M[FWD])
# AIC
# calculate the long way for M1
ll1 <- dnorm(strike$Strike, mean = predict(M1),
             sd = summary(M1)$sigma*sqrt(M1$df/nrow(strike)), # mle divides by n
             log = TRUE) # loglikelihood
ll1 <- sum(ll1)
p1 <- length(M1$coef)+1 # number of parameters
AIC1 <- -2*ll1 + 2*p1
# calculate the R way
AIC1 - AIC(M1)
AIC2 <- AIC(M2) # for M2
c(AIC1 = AIC1, AIC2 = AIC2) # favors Mfwd
# PRESS
press1 <- resid(M1)/(1-hatvalues(M1)) # M1
press2 <-  resid(M2)/(1-hatvalues(M2)) # M2
c(PRESS1 = sum(press1^2), PRESS2 = sum(press2^2)) # favors Mfwd

# Cross-validation
nreps <- 2e3 # number of replications
ntot <- nrow(strike) # total number of observations
ntrain <- 500 # size of training set
ntest <- ntot-ntrain # size of test set
sse1 <- rep(NA, nreps) # sum-of-square errors for each CV replication
sse2 <- rep(NA, nreps)
ssefwd <- rep(NA, nreps)
Lambda <- rep(NA, nreps) # likelihod ratio statistic for each replication
system.time({
  for(ii in 1:nreps) {
    if(ii%%400 == 0) message("ii = ", ii)
    # randomly select training observations
    train.ind <- sample(ntot, ntrain) # training observations
    # this is the more straightforward way of calculating the CV parameter estimates
    # M1.cv <- lm(Strike ~ read + prog + race + ses + locus + read:prog + prog:ses,
    #             data = strike, subset = train.ind)
    # this is the faster R way
    M1.cv <- update(M1, subset = train.ind)
    M2.cv <- update(M2, subset = train.ind)
    Mfwd.cv <- update(Mfwd, subset = train.ind)
    # testing residuals for both models
    # that is, testing data - predictions with training parameters
    M1.res <- strike$Strike[-train.ind] - predict(M1.cv, newdata = strike[-train.ind,])
    M2.res <- strike$Strike[-train.ind] - predict(M2.cv, newdata = strike[-train.ind,])
    Mfwd.res <- strike$Strike[-train.ind] - predict(Mfwd.cv, newdata = strike[-train.ind,])
    # total sum of square errors
    sse1[ii] <- sum((M1.res)^2)
    sse2[ii] <- sum((M2.res)^2)
    ssefwd[ii] <- sum((Mfwd.res)^2)
    # testing likelihood ratio
    M1.sigma <- sqrt(sum(resid(M1.cv)^2)/ntrain) # MLE of sigma
    M2.sigma <- sqrt(sum(resid(M2.cv)^2)/ntrain)
    Lambda[ii] <- sum(dnorm(M1.res, mean = 0, sd = M1.sigma, log = TRUE))
    Lambda[ii] <- Lambda[ii] - sum(dnorm(M2.res, mean = 0, sd = M2.sigma, log = TRUE))
  } })

message("1:M1, 2:M2, 3:Mfwd")
c(SSE1 = mean(sse1), SSE2 = mean(sse2), SSEFWD = mean(ssefwd)) # favors Mfwd
mean(Lambda) # log(Lambda) < 0 which favors Mfwd