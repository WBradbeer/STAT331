# load file
strike <- read.csv("C:/Users/Karthiga/Desktop/STAT 331/project/strikes_clean.csv")
# encode categoricals as Factors
strike$Centr <- as.factor(strike$Centr)
strike$Country <- as.factor(strike$Country)

# Investigate to look for patterns
summary(strike)
pairs(strike)
#correlation matrix
d <- data.frame(strike$Strike,
                strike$Unemp,
                strike$Infl,
                strike$Demo,
                strike$Dens)
cor(d)

# Reason that these should affect strike
M1 = lm(Strike ~ Infl*Centr*Unemp, data = strike)
# anova shows unemp doesn't add much 
anova(M1)

# remove Unemp from M1 to get current model that works well 
M2 = lm(Strike ~ Infl*Centr, data = strike)
# basic model
M0 <- lm(Strike ~ 1, data = strike)
# full model without country - interaction effects
Mfull <- lm(Strike ~ (. - Country)^2 + Country, data = strike)
# forward step model from intercept to full
Mfwd <- step(M0, scope=list(lower=M0, upper=Mfull), direction='forward', trace=FALSE)
# backward step model from full to intercept
Mback <- step(Mfull, scope=list(lower=M0, upper=Mfull), direction='backward', trace=FALSE)
# step model starting from model grenerated from exploration
Mstep <- step(M2, scope=list(lower=M0, upper=Mfull), direction='both', trace=FALSE)

models <- list(M0, M1, M2, Mfwd, Mback, Mstep)

n <- length(strike$Strike)
samp.size <- 400
results <- rep(0, length(models))

for(i in 1:100) {
  
  # seperate training and test data
  samp <- sample.int(n, size = samp.size)
  strike.train <- strike[samp,]
  strike.test <- strike[-samp,]
  j <- 1
  for(m in models){
    # train model
    m.trn <- update(m, data=strike.train)
    # test model
    results[j] <- results[j] + sum((predict.lm(m, strike.test) - strike.test$Strike)^2)
    j <- j+1
  }
}


#compare 3 models
Mfwd$call
Mback$call
Mstep$call

Mnames <- expression(M[FWD], M[BACK])

# AIC
AIC1 <- AIC(Mfwd) # for Mfwd
AIC2 <- AIC(Mback) # for Mback
c(AIC1 = AIC1, AIC2 = AIC2) #favours Mback

# PRESS
press1 <- resid(Mfwd)/(1-hatvalues(Mfwd)) # M1
press2 <-  resid(Mback)/(1-hatvalues(Mback)) # M2
c(PRESS1 = sum(press1^2), PRESS2 = sum(press2^2)) # favors M2

#Cross-validation
nreps <- 2e3 # number of replications
ntot <- nrow(strike) # total number of observations
ntrain <- 500 # size of training set
ntest <- ntot-ntrain # size of test set
sse1 <- rep(NA, nreps) # sum-of-square errors for each CV replication
sse2 <- rep(NA, nreps)
Lambda <- rep(NA, nreps) # likelihod ratio statistic for each replication

system.time({
  for(ii in 1:nreps) {
    if(ii%%400 == 0) message("ii = ", ii)
    # randomly select training observations
    train.ind <- sample(ntot, ntrain) # training observations
    # this is the more straightforward way of calculating the CV parameter estimates
    # Mfwd.cv <- lm(math ~ read + prog + race + ses + locus + read:prog + prog:ses,
    #             data = strike, subset = train.ind)
    # this is the faster R way
    Mfwd.cv <- update(Mfwd, subset = train.ind)
    Mback.cv <- update(Mback, subset = train.ind)
    # testing residuals for both models
    # that is, testing data - predictions with training parameter
    Mfwd.res <- strike$math[-train.ind] - predict(Mfwd.cv, newdata = strike[-train.ind,])
    Mback.res <- strike$math[-train.ind] - predict(Mback.cv, newdata = strike[-train.ind,])
    # total sum of square errors
    sse1[ii] <- sum((Mfwd.res)^2)
    sse2[ii] <- sum((Mback.res)^2)
    # testing likelihood ratio
    Mfwd.sigma <- sqrt(sum(resid(Mfwd.cv)^2)/ntrain) # MLE of sigma
    Mback.sigma <- sqrt(sum(resid(Mback.cv)^2)/ntrain)
    Lambda[ii] <- sum(dnorm(Mfwd.res, mean = 0, sd = Mfwd.sigma, log = TRUE))
    Lambda[ii] <- Lambda[ii] - sum(dnorm(Mback.res, mean = 0, sd = Mback.sigma, log = TRUE))
  }
})



c(SSE1 = mean(sse1), SSE2 = mean(sse2)) # favors M2
mean(Lambda) # log(Lambda) < 0 which favors M2

#plot PRESS statistics
par(mar = c(3, 6, 1, 1))
boxplot(x = list(press1^2, press2^2), names = Mnames,
        ylab = expression(PRESS[i]^2), col = c("yellow", "orange"))

#plot cross-validation SSE and Lambda
par(mfrow = c(1,2))
par(mar = c(4.5, 4.5, .1, .1))
boxplot(x = list(sse1, sse2), names = Mnames, cex = .7,
        ylab = expression(SS[err]^{test}), col = c("yellow", "orange"))
hist(Lambda, breaks = 50, freq = FALSE, xlab = expression(Lambda^{test}),
     main = "", cex = .7)
abline(v = mean(Lambda), col = "red") # average value


#plot(strike$Dens, log(strike$Strike), col = strike$Country)
#plot(log(strike$Unemp), log(strike$Strike), col = strike$Country)
#legend(93,13, unique(strike$Country),col=1:length(strike$Country),pch=1)