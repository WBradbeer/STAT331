
# load file
strike <- read.csv("strikes_clean.csv")
# encode categoricals as Factors
strike$Centr <- as.factor(strike$Centr)
strike$Country <- as.factor(strike$Country)
#strike$Unemp2 <- strike$Unemp^2
#strike$Demo2 <- strike$Demo^2
#strike$Infl2 <- strike$Infl^2

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


models <- list(M2, )

 
n <- length(strike$Strike)
samp.size <- 400
results <- rep(0, length(models))
results2 <- rep(0, length(models))

for(i in 1:2000) {
  # seperate training and test data
  train.ind <- sample(n, samp.size)
  strike.train <- strike[train.ind,]
  strike.test <- strike[-train.ind,]
  
  j <- 1
  for(m in models){
    # train model
    m.trn <- update(m, subset=train.ind)
    # test model
    m.res = predict.lm(m.trn, strike.test) - strike.test$Strike
    results[j] <- results[j] + sum(m.res^2)
    j <- j+1
  }
}

# model diagnostics
res <- resid(Mfwd)
H <-  model.matrix(Mfwd)
H <- H %*% solve(crossprod(H), t(H))
h <- diag(H)
# studentized residuals
res.stu <- res/sqrt(1-h)

# standardized student residuals
res.stand <- res.stu/sigma(Mfwd)

# diagnostic plots
# scatter plot of residuals vs predictions
cex <- 0.8 
par(mar = c(4,4,0.1,0.1))
plot(predict(Mfwd), res, pch=21, bg="black", cex=cex, cex.axis=cex, 
     xlab="Predicted # of strike-hours", ylab="Residual strike-hours")
points(predict(Mfwd), res.stu, pch=21, bg="red", cex=cex)
legend(x="topleft", c("Residuals", "Studentized Residuals"), pch=21,
       pt.bg = c("black", "red"), pt.cex=cex, cex=cex)
# seems like multiplicative error

# histogram of studentized residuals
par(mar = c(4,4,0.1,0.1))
hist.plot <- hist(res.stand, breaks=100, freq = FALSE, cex.axis=cex, main="",
                  xlab="Standardized, studentized residual strike activity")
curve(dnorm(x), col='red', add=TRUE)

# leverage


