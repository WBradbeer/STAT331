# load file
strike <- read.csv("strikes_clean.csv")
# encode categoricals as Factors
strike$Centr <- as.factor(strike$Centr)
strike$Country <- as.factor(strike$Country)
#strike$Unemp2 <- strike$Unemp^2
#strike$Demo2 <- strike$Demo^2
strike$Infl2 <- strike$Infl^2

# Investigate to look for patterns  
summary(strike)
pairs(strike)


# Reason that these should affect strike
M1 = lm(Strike ~ Infl*Centr*Unemp, data = strike)
# anova shows unemp doesn't add much 
anova(M1)
# remove Unemp from M1 to get current model that works well 
M2 = lm(Strike ~ Infl*Centr, data = strike)

M3 = lm(Strike ~ Infl2*Centr, data = strike)

M4 = lm(Strike ~ Country*Year, data = strike)

M5 = lm(Strike ~ Infl*Centr + Country, data = strike)

M6 = lm(Strike ~ Centr + Infl + log(Dens) + log(Demo) + Centr:Infl + Centr:log(Dens) + Infl:log(Demo), data = strike)


# basic model
M0 <- lm(Strike ~ 1, data = strike)
# full model without country - interaction effects
Mfull <- lm(Strike ~ (. - Country - Year)^2 + Country, data = strike)
# forward step model from intercept to full
Mfwd <- step(M0, scope=list(lower=M0, upper=Mfull), direction='forward', trace=FALSE)
# backward step model from full to intercept
Mback <- step(Mfull, scope=list(lower=M0, upper=Mfull), direction='backward', trace=FALSE)
# step model starting from model grenerated from exploration
Mstep <- step(M2, scope=list(lower=M0, upper=Mfull), direction='both', trace=FALSE)

M1 <- M2
M2 <- Mfwd
Mnames <- expression(M[2], M[fwd])
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
    # testing residuals for both models
    # that is, testing data - predictions with training parameters
    M1.res <- strike$Strike[-train.ind] - predict(M1.cv, newdata = strike[-train.ind,])
    M2.res <- strike$Strike[-train.ind] - predict(M2.cv, newdata = strike[-train.ind,])
    # total sum of square errors
    sse1[ii] <- sum((M1.res)^2)
    sse2[ii] <- sum((M2.res)^2)
    # testing likelihood ratio
    M1.sigma <- sqrt(sum(resid(M1.cv)^2)/ntrain) # MLE of sigma
    M2.sigma <- sqrt(sum(resid(M2.cv)^2)/ntrain)
    Lambda[ii] <- sum(dnorm(M1.res, mean = 0, sd = M1.sigma, log = TRUE))
    Lambda[ii] <- Lambda[ii] - sum(dnorm(M2.res, mean = 0, sd = M2.sigma, log = TRUE))
  } })

message("1:M1, 2:Mfwd")
c(SSE1 = mean(sse1), SSE2 = mean(sse2)) # favors smaller of 2
# in units of strike days
c(SSE1 = sqrt(mean(sse1)/ntest), SSE2 = sqrt(mean(sse2)/ntest)) # favors smaller of 2
mean(Lambda) # log(Lambda) < 0 which favors M2 > 0 M1


dev.off()

# model diagnostics
res1 <- resid(M1)
H1 <-  model.matrix(M1)
H1 <- H1 %*% solve(crossprod(H1), t(H1))
h1 <- diag(H1)
# studentized residuals
res.stu1 <- res1/sqrt(1-h1)

# standardized student residuals
res.stand1 <- res.stu1/sigma(M1)

cex <- 0.8 
par(mar = c(4,4,0.1,0.1))
plot(predict(M1), res1, pch=21, bg="black", cex=cex, cex.axis=cex,
     xlab="Predicted # of strike-hours", ylab="Residual strike-hours")
points(predict(M1), res.stu1, pch=21, bg="red", cex=cex)
legend(x="topleft", c("Residuals", "Studentized Residuals"), pch=21,
       pt.bg = c("black", "red"), pt.cex=cex, cex=cex)


res2 <- resid(M2)
H2 <- model.matrix(M2)
H2 <- H2 %*% solve(crossprod(H2), t(H2))
h2 <- diag(H2)
# studentized residuals
res.stu2 <- res2/sqrt(1-h2)

# standardized student residuals
res.stand2 <- res.stu2/sigma(M2)

cex <- 0.8 
par(mar = c(4,4,0.1,0.1))
plot(predict(M2), res, pch=21, bg="black", cex=cex, cex.axis=cex,  
     xlab="Predicted # of strike-hours", ylab="Residual strike-hours")
points(predict(M2), res.stu, pch=21, bg="red", cex=cex)
legend(x="topleft", c("Residuals", "Studentized Residuals"), pch=21,
       pt.bg = c("black", "red"), pt.cex=cex, cex=cex)
# seems like multiplicative error

# histogram of studentized residuals
par(mar = c(4,4,0.1,0.1))
hist.plot1 <- hist(res.stand1, breaks=75, freq = FALSE, cex.axis=cex, main="",
                  xlab="Standardized, studentized residual strike activity")
curve(dnorm(x), col='red', add=TRUE)

# histogram of studentized residuals
par(mar = c(4,4,0.1,0.1))
hist.plot2 <- hist(res.stand2, breaks=100, freq = FALSE, cex.axis=cex, main = "",
                  xlab="Standardized, studentized residual strike activity")
curve(dnorm(x), col='red', add=TRUE)

#plot cross-validation SSE and Lambda
par(mfrow = c(1,2))
par(mar = c(4.5, 4.5, .1, .1))
boxplot(x = list(sqrt(sse1/ntest), sqrt(sse2/ntest)), names = Mnames, cex = .7,
        ylab = expression(RMSE^{test}), col = c("yellow", "orange"))
hist(Lambda, breaks = 50, freq = FALSE, xlab = expression(Lambda^{test}),
     main = "", cex = .7)
abline(v = mean(Lambda), col = "red") # average value
