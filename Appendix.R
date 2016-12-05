# load file
strike <- read.csv("Documents/Waterloo/4A/331/project/strikes_clean.csv")
# encode categoricals as Factors
strike$Centr <- as.factor(strike$Centr)
strike$Country <- as.factor(strike$Country)


# Section 2: Model Selection


# Investigate to look for patterns  
summary(strike)
pairs(strike)


# Reason that these should affect strike
M1 <- lm(Strike ~ Infl*Centr*Unemp, data = strike)
# anova shows Unemp doesn't add much 
anova(M1)
# Remove Unemp from M1 to get current model that works well
M2 <- lm(Strike ~ Infl*Centr, data = strike)


# Section 2.1: Automated Model Selection


# M0: Basic model
M0 <- lm(Strike ~ 1, data = strike)
# Mfull: Full model without Centr and without Country interaction effects
Mfull <- lm(Strike ~ (. - Country - Centr)^2 + Country, data = strike)
# Remove Centr from M1 since Mfull does not have Centr
M3 <- lm(Strike ~ Infl*Unemp, data = strike)


# Mfwd: Forward step model from intercept to full
Mfwd <- step(M0, scope=list(lower=M0, upper=Mfull), direction='forward', trace=FALSE)
# Mback: Backward step model from full to intercept
Mback <- step(Mfull, scope=list(lower=M0, upper=Mfull), direction='backward', trace=FALSE)
# Mstep: Step model starting from M3 from above
Mstep <- step(M3, scope=list(lower=M0, upper=Mfull), direction='both', trace=FALSE)


# Run summary() to see call, residuals and coefficients of three automated models
summary(Mfwd)
summary(Mback)
summary(Mstep)








# Section 2.2: Automated Model Selection with Macroeconomic Variables


# M0 and M2 remain same as above
# Mfullmacro: Full model without Country or Year
Mfullmacro <- lm(Strike ~ (. - Country - Year)^2, data = strike)


# Mfwdmacro: Forward step model from intercept to full
Mfwdmacro <- step(M0, scope=list(lower=M0, upper=Mfullmacro), direction='forward', trace=FALSE)
# Mbackmacro: Backward step model from full to intercept
Mbackmacro <- step(Mfullmacro, scope=list(lower=M0, upper=Mfullmacro), direction='backward', trace=FALSE)
# Mstepmacro: Step model starting from M2 from above
Mstepmacro <- step(M2, scope=list(lower=M0, upper=Mfullmacro), direction='both', trace=FALSE)


# Run summary() to see call, residuals and coefficients of three automated models
summary(Mfwdmacro)
summary(Mbackmacro)
summary(Mstepmacro)


# Section 2.2: Automated Model Selection with Macroeconomic Variables


# Use Centr as a numeric to test correlation
centrNum <- as.numeric(levels(strike$Centr))[strike$Centr]
# Load all macroeconomic variables into correlation matrix, exclude Country and Year
correlMatrix <- cbind(strike$Unemp, strike$Infl, strike$Demo, strike$Dens, centrNum)
# Output correlation between each variable in correlation matrix
cor(correlMatrix)


# Update correlation matrix to have log of Demo and Dens
correlMatrix <- cbind(strike$Unemp, strike$Infl, log(strike$Demo), log(strike$Dens), centrNum)
# Output correlation between each variable in correlation matrix
cor(correlMatrix)


# Automated Model Selection: M0 and M2 remain same as above
# Mfulllog: Full model without Country or Year and with log of Demo and Dens
Mfulllog <- lm(Strike ~ (. - Country - Year - Demo - Dens + log(Demo) + log(Dens))^2, data = strike)


# Mlogmacro: Step model starting from M2 from above
Mlogmacro <- step(M2, scope=list(lower=M0, upper=Mfulllog), direction='both', trace=FALSE)


# Run summary() to see call, residuals and coefficients of three automated models
summary(Mlogmacro)









# Section 3: Model Diagnostics

# For model diagnostics we will be calling M2 from above M_a
# and we will be calling Mfwd from above M_b 
M_a <- M2 
M_b <- Mfwd
Mnames <- expression(M[2], M[fwd])

# Section 3.1: Residual Plots

# First we find the residual plot of M_a
res_a <- resid(M_a)
H_a <-  model.matrix(M_a)
H_a <- H_a %*% solve(crossprod(H_a), t(H_a))
h_a <- diag(H_a)
# Studentized residuals
res.stu_a <- res_a/sqrt(1-h_a)

# Standardized student residuals
res.sta_a = res_a/sigma(M_a)
res.stastu_a <- (res_a/sigma(M_a))/sqrt(1-h_a)

cex <- 0.8 
par(mar = c(4,4,0.1,0.1))
plot(predict(M_a), res_a, pch=21, bg="black", cex=cex, cex.axis=cex,
     xlab="Predicted # of strike-hours", ylab="Residual strike-hours")
points(predict(M_a), res.stu_a, pch=21, bg="red", cex=cex)
legend(x="topleft", c("Residuals", "Studentized Residuals"), pch=21,
       pt.bg = c("black", "red"), pt.cex=cex, cex=cex)


# Next we find the residual plot of M_b
res_b <- resid(M_b)
H_b <- model.matrix(M_b)
H_b <- H_b %*% solve(crossprod(H_b), t(H_b))
h_b <- diag(H_b)
# Studentized residuals
res.stu_b <- res_b/sqrt(1-h_b)

# Standardized student residuals
res.sta_b = res_b/sigma(M_b)
res.stastu_b <- (res_b/sigma(M_b))/sqrt(1-h_b)

cex <- 0.8 
par(mar = c(4,4,0.1,0.1))
plot(predict(M_b), res_b, pch=21, bg="black", cex=cex, cex.axis=cex,  
     xlab="Predicted # of strike-hours", ylab="Residual strike-hours")
points(predict(M_b), res.stu_b, pch=21, bg="red", cex=cex)
legend(x="topleft", c("Residuals", "Studentized Residuals"), pch=21,
       pt.bg = c("black", "red"), pt.cex=cex, cex=cex)



# Next, we will plot the histogram of studentized residuals for M_a
par(mar = c(4,4,0.1,0.1))
hist.plot_a <- hist(res.stastu_a, breaks=75, freq = FALSE, cex.axis=cex, main="",
                   xlab="Standardized, studentized residual strike activity")
curve(dnorm(x), col='red', add=TRUE)

# Next, we will plot the histogram of studentized residuals for M_b
par(mar = c(4,4,0.1,0.1))
hist.plot_b <- hist(res.stastu_b, breaks=100, freq = FALSE, cex.axis=cex, main = "",
                   xlab="Standardized, studentized residual strike activity")
curve(dnorm(x), col='red', add=TRUE)
















# Section 3.2: Leverage and Influence Measures


#Leverage and Influence Measures of M_a

head(diag(H_a)) 
hat_a <- hatvalues(M_a) # the R way
range(hat_a - diag(H_a)) 
y.hat_a <- predict(M_a)
sigma.hat_a <- summary(M_a)$sigma

# PRESS residuals
press_a <- res_a/(1-hat_a)
# DFFITS residuals
dfts_a <- dffits(M_a) # the R way 

# standardize each of these suchthat they are identical at the average leverage value
p_a <- length(coef(M_a))
n_a <- nobs(M_a)
hbar_a <- p_a/n_a # average leverage
res.stastu_a <- res.stastu_a*sqrt(1-hbar_a) # at hat_a = hbar_a, res.stastu_a = res.sta_a
press_a <- press_a*(1-hbar_a)/sigma.hat_a # at hat_a = hbar_a, press_a = res.sta_a
dfts_a <- dfts_a*(1-hbar_a)/sqrt(hbar_a)# at hat_a = hbar_a, dfts_a = res.sta_a 

clrs <- rep("black", len = n_a)
strike_a = data.frame(strike$Strike, strike$Centr, strike$Infl)
pairs(strike_a, pch = 21, bg = clrs)



# plot all residuals
par(mfrow = c(1,2), mar = c(4,4,.1,.1)) # against predicted values
cex <- .8
plot(y.hat_a, rep(0, length(y.hat_a)), type = "n", # empty plot to get the axis range
     ylim = range(res.sta_a, res.stastu_a, press_a, dfts_a), cex.axis = cex,
     xlab = "Predicted Values", ylab = "Residuals")
# dotted line connecting each observations residuals for better visibility 
segments(x0 = y.hat_a,
         y0 = pmin(res.sta_a, res.stastu_a, press_a, dfts_a),
         y1 = pmax(res.sta_a, res.stastu_a, press_a, dfts_a),
         lty = 2)
points(y.hat_a, res.sta_a, pch = 21, bg = "black", cex = cex)
points(y.hat_a, res.stastu_a, pch = 21, bg = "blue", cex = cex)
points(y.hat_a, press_a, pch = 21, bg = "red", cex = cex)
points(y.hat_a, dfts_a, pch = 21, bg= "orange", cex = cex)
# against leverages
plot(hat_a, rep(0, length(y.hat_a)), type = "n", cex.axis = cex,
     ylim = range(res.sta_a, res.stastu_a, press_a, dfts_a),
     xlab = "Leverages", ylab = "Residuals")
segments(x0 = hat_a,
         y0 = pmin(res.sta_a, res.stastu_a, press_a, dfts_a),
         y1 = pmax(res.sta_a, res.stastu_a, press_a, dfts_a),
         lty = 2)
points(hat_a, res.sta_a, pch = 21, bg= "black", cex = cex)
points(hat_a, res.stastu_a, pch = 21, bg= "blue", cex = cex)
points(hat_a, press_a, pch = 21, bg = "red", cex = cex)
points(hat_a, dfts_a, pch = 21, bg = "orange", cex = cex)
abline(v = hbar_a, col = "grey60",lty = 2)
legend("topright", legend = c("Standardized", "Studentized", "PRESS", "DFFITS"),
       pch = 21, pt.bg = c("black", "blue", "red", "orange"), title = "Residual Type:",
       cex = cex, pt.cex = cex,)

# cook's distance vs. leverage
D_a <- cooks.distance(M_a)
# flag some of the points
infl.ind <- which.max(D_a) # top influence point

lev.ind <- hat_a > 2*hbar_a # leveragemore than 2x the average
clrs[lev.ind] <- "blue"
clrs[infl.ind] <- "red"
par(mfrow = c(1,1), mar = c(4,4,1,1))
cex <- .8
plot(hat_a, D_a, xlab = "Leverage", ylab = "Cook's Influence Measure",
     pch = 21, bg = clrs, cex = cex, cex.axis = cex)
p_a <- length(coef(M_a))
n_a <- nrow(strike)
hbar_a <- p_a/n_a # average leverage
abline(v = 2*hbar_a, col = "grey60", lty = 2) # 2x average leverage
legend("topleft", legend = c("High Leverage", "High Influence"), pch = 21,
       pt.bg = c("blue", "red"),cex = cex, pt.cex = cex)

# pairs plot of all covariates to see high leverage & influence
pairs(strike_a, pch = 21, bg = clrs)

# predictions with omitted values
omit.ind <- c(infl.ind, # most influential
              which.max(hat_a)) # highest leverage
names(omit.ind) <- c("infl", "lev")
omit.ind




#Leverage and Influence Measures of M_b

head(diag(H_b)) 
hat_b <- hatvalues(M_b) # the R way
range(hat_b - diag(H_b)) 
y.hat_b <- predict(M_b)
sigma.hat_b <- summary(M_b)$sigma

# PRESS residuals
press_b <- res_b/(1-hat_b)
# DFFITS residuals
dfts_b <- dffits(M_b) # the R way 

# standardize each of these suchthat they are identical at the average leverage value
p_b <- length(coef(M_b))
n_b <- nobs(M_b)
hbar_b <- p_b/n_b # average leverage
res.stastu_b <- res.stastu_b*sqrt(1-hbar_b) # at hat_b = hbar_b, res.stastu_b = res.sta_b
press_b <- press_b*(1-hbar_b)/sigma.hat_b # at hat_b = hbar_b, press_b = res.sta_b
dfts_b <- dfts_b*(1-hbar_b)/sqrt(hbar_b)# at hat_b = hbar_b, dfts_b = res.sta_b 

clrs <- rep("black", len = n_b)
strike_b = data.frame(strike$Strike, strike$Country, strike$Infl, strike$Unemp)
pairs(strike_b, pch = 21, bg = clrs)



# plot all residuals
par(mfrow = c(1,2), mar = c(4,4,.1,.1)) # against predicted values
cex <- .8
plot(y.hat_b, rep(0, length(y.hat_b)), type = "n", # empty plot to get the axis range
     ylim = range(res.sta_b, res.stastu_b, press_b, dfts_b), cex.axis = cex,
     xlab = "Predicted Values", ylab = "Residuals")
# dotted line connecting each observations residuals for better visibility 
segments(x0 = y.hat_b,
         y0 = pmin(res.sta_b, res.stastu_b, press_b, dfts_b),
         y1 = pmax(res.sta_b, res.stastu_b, press_b, dfts_b),
         lty = 2)
points(y.hat_b, res.sta_b, pch = 21, bg = "black", cex = cex)
points(y.hat_b, res.stastu_b, pch = 21, bg = "blue", cex = cex)
points(y.hat_b, press_b, pch = 21, bg = "red", cex = cex)
points(y.hat_b, dfts_b, pch = 21, bg= "orange", cex = cex)
# against leverages
plot(hat_b, rep(0, length(y.hat_b)), type = "n", cex.axis = cex,
     ylim = range(res.sta_b, res.stastu_b, press_b, dfts_b),
     xlab = "Leverages", ylab = "Residuals")
segments(x0 = hat_b,
         y0 = pmin(res.sta_b, res.stastu_b, press_b, dfts_b),
         y1 = pmax(res.sta_b, res.stastu_b, press_b, dfts_b),
         lty = 2)
points(hat_b, res.sta_b, pch = 21, bg= "black", cex = cex)
points(hat_b, res.stastu_b, pch = 21, bg= "blue", cex = cex)
points(hat_b, press_b, pch = 21, bg = "red", cex = cex)
points(hat_b, dfts_b, pch = 21, bg = "orange", cex = cex)
abline(v = hbar_b, col = "grey60",lty = 2)
legend("topright", legend = c("Standardized", "Studentized", "PRESS", "DFFITS"),
       pch = 21, pt.bg = c("black", "blue", "red", "orange"), title = "Residual Type:",
       cex = cex, pt.cex = cex,)

# cook's distance vs. leverage
D_b <- cooks.distance(M_b)
# flag some of the points
infl.ind <- which.max(D_b) # top influence point

lev.ind <- hat_b > 2*hbar_b # leveragemore than 2x the average
clrs[lev.ind] <- "blue"
clrs[infl.ind] <- "red"
par(mfrow = c(1,1), mar = c(4,4,1,1))
cex <- .8
plot(hat_b, D_b, xlab = "Leverage", ylab = "Cook's Influence Measure",
     pch = 21, bg = clrs, cex = cex, cex.axis = cex)
p_b <- length(coef(M_b))
n_b <- nrow(strike)
hbar_b <- p_b/n_b # average leverage
abline(v = 2*hbar_b, col = "grey60", lty = 2) # 2x average leverage
legend("topleft", legend = c("High Leverage", "High Influence"), pch = 21,
       pt.bg = c("blue", "red"),cex = cex, pt.cex = cex)

# pairs plot of all covariates to see high leverage & influence
pairs(strike_b, pch = 21, bg = clrs)

# predictions with omitted values
omit.ind <- c(infl.ind, # most influential
              which.max(hat_b)) # highest leverage
names(omit.ind) <- c("infl", "lev")
omit.ind



press_a <-  resid(M_a)/(1-hatvalues(M_a)) # M2
press_b <- resid(M_b)/(1-hatvalues(M_b)) # Mfwd
c(PRESS_a = sum(press_a^2), PRESS_b = sum(press_b^2)) # favors M2




# Section 3.3: Cross-Validation

# Cross-validation
nreps <- 2e3 # number of replications
ntot <- nrow(strike) # total number of observations
ntrain <- 500 # size of training set
ntest <- ntot-ntrain # size of test set
sse_a <- rep(NA, nreps) # sum-of-square errors for each CV replication
sse_b <- rep(NA, nreps)
Lambda <- rep(NA, nreps) # likelihod ratio statistic for each replication
system.time({
  for(ii in 1:nreps) {
    if(ii%%400 == 0) message("ii = ", ii)
    # randomly select training observations
    train.ind <- sample(ntot, ntrain) # training observations
    # this is the more straightforward way of calculating the CV parameter estimates
    # M_a.cv <- lm(Strike ~ read + prog + race + ses + locus + read:prog + prog:ses,
    #             data = strike, subset = train.ind)
    # this is the faster R way
    M_a.cv <- update(M_a, subset = train.ind)
    M_b.cv <- update(M_b, subset = train.ind)
    # testing residuals for both models
    # that is, testing data - predictions with training parameters
    M_a.res <- strike$Strike[-train.ind] - predict(M_a.cv, newdata = strike[-train.ind,])
    M_b.res <- strike$Strike[-train.ind] - predict(M_b.cv, newdata = strike[-train.ind,])
    # total sum of square errors
    sse_a[ii] <- sum((M_a.res)^2)
    sse_b[ii] <- sum((M_b.res)^2)
    # testing likelihood ratio
    M_a.sigma <- sqrt(sum(resid(M_a.cv)^2)/ntrain) # MLE of sigma
    M_b.sigma <- sqrt(sum(resid(M_b.cv)^2)/ntrain)
    Lambda[ii] <- sum(dnorm(M_a.res, mean = 0, sd = M_a.sigma, log = TRUE))
    Lambda[ii] <- Lambda[ii] - sum(dnorm(M_b.res, mean = 0, sd = M_b.sigma, log = TRUE))
  } })

message("1:M2, 2:Mfwd")
c(SSE_a = mean(sse_a), SSE_b = mean(sse_b)) # favors smaller of 2
# in units of strike days
c(SSE_a = sqrt(mean(sse_a)/ntest), SSE_b = sqrt(mean(sse_b)/ntest)) # favors smaller of 2
mean(Lambda) # log(Lambda)


#plot cross-validation SSE and Lambda
par(mfrow = c(1,2))
par(mar = c(4.5, 4.5, .1, .1))
boxplot(x = list(sqrt(sse_a/ntest), sqrt(sse_b/ntest)), names = Mnames, cex = .7,
        ylab = expression(RMSE^{test}), col = c("yellow", "orange"))
hist(Lambda, breaks = 50, freq = FALSE, xlab = expression(Lambda^{test}),
     main = "", cex = .7)
abline(v = mean(Lambda), col = "red") # average value