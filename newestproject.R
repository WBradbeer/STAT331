# load file
strike <- read.csv("C:/Users/Karthiga/Desktop/STAT 331/project/strikes_clean.csv")
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
# res = original residuals
# res.stand = standardized student residuals
# res.stu = studentized residuals

# H = hat matrix
# h = diagonal of hat matrix H

hval <- hatvalues(Mfwd)
range(hval - h) 

# studentized residuals
#res.stu

# PRESS residuals
press <- res/(1-h)

# DFFITS residuals
dfts <- dffits(M)



#Simple macro
M2 <- lm(formula = Strike ~ Infl * Centr, data = strike)
#Country Fwd model
Mfwd <- lm(formula = Strike ~ Country + Infl + Unemp, data = strike)


M <- M2
Msum <- summary(M) 
# residuals vs predicted values
y.hat <- predict(M) # predicted values
sigma.hat <- Msum$sigma
res <- resid(M)
# original residuals$
stan.res <- res/sigma.hat # standardized residuals # compute leverages
X <- model.matrix(M)
H <- X %*% solve(crossprod(X), t(X)) # HAT matrix
head(diag(H)) 

h <- hatvalues(M) # the R way
range(h - diag(H)) 

# studentized residuals
stud.res <- stan.res/sqrt(1-h)
# PRESS residuals
press <- res/(1-h)
# DFFITS residuals
dfts <- dffits(M) # the R way 

# standardize each of these suchthat they are identical at the average leverage value
p <- length(coef(M))
n <- nobs(M)
hbar <- p/n # average leverage
stud.res <- stud.res*sqrt(1-hbar) # at h = hbar, stud.res = stan.res
press <- press*(1-hbar)/sigma.hat # at h = hbar, press = stan.res
dfts <- dfts*(1-hbar)/sqrt(hbar)# at h = hbar, dfts = stan.res 

# plot all residuals
par(mfrow = c(1,2), mar = c(4,4,.1,.1)) # against predicted values
cex <- .8
plot(y.hat, rep(0, length(y.hat)), type = "n", # empty plot to get the axis range
     ylim = range(stan.res, stud.res, press, dfts), cex.axis = cex,
     xlab = "Predicted Values", ylab = "Residuals")
# dotted line connecting each observations residuals for better visibility 
segments(x0 = y.hat,
         y0 = pmin(stan.res, stud.res, press, dfts),
         y1 = pmax(stan.res, stud.res, press, dfts),
         lty = 2)
points(y.hat, stan.res, pch = 21, bg = "black", cex = cex)
points(y.hat, stud.res, pch = 21, bg = "blue", cex = cex)
points(y.hat, press, pch = 21, bg = "red", cex = cex)
points(y.hat, dfts, pch = 21, bg= "orange", cex = cex)
# against leverages
plot(h, rep(0, length(y.hat)), type = "n", cex.axis = cex,
     ylim = range(stan.res, stud.res, press, dfts),
     xlab = "Leverages", ylab = "Residuals")
segments(x0 = h,
         y0 = pmin(stan.res, stud.res, press, dfts),
         y1 = pmax(stan.res, stud.res, press, dfts),
         lty = 2)
points(h, stan.res, pch = 21, bg= "black", cex = cex)
points(h, stud.res, pch = 21, bg= "blue", cex = cex)
points(h, press, pch = 21, bg = "red", cex = cex)
points(h, dfts, pch = 21, bg = "orange", cex = cex)
abline(v = hbar, col = "grey60",lty = 2)
legend("topright", legend = c("Standardized", "Studentized", "PRESS", "DFFITS"),
       pch = 21, pt.bg = c("black", "blue", "red", "orange"), title = "Residual Type:",
       cex = cex, pt.cex = cex,)

# cook's distance vs. leverage
D <- cooks.distance(M)
# flag some of the points
infl.ind <- which.max(D) # top influence point

lev.ind <- h > 2*hbar # leveragemore than 2x the average
clrs <- rep("black", len = n)
clrs[lev.ind] <- "blue"
clrs[infl.ind] <- "red"
par(mfrow = c(1,1), mar = c(4,4,1,1))
cex <- .8
plot(h, D, xlab = "Leverage", ylab = "Cook's Influence Measure",
     pch = 21, bg = clrs, cex = cex, cex.axis = cex)
p <- length(coef(M))
n <- nrow(strike)
hbar <- p/n # average leverage
abline(v = 2*hbar, col = "grey60", lty = 2) # 2x average leverage
legend("topleft", legend = c("High Leverage", "High Influence"), pch = 21,
       pt.bg = c("blue", "red"),cex = cex, pt.cex = cex)

# pairs plot of all covariates to see high leverage & influence
pairs(strikeMfwd, pch = 21, bg = clrs)

# predictions with omitted values
omit.ind <- c(infl.ind, # most influential
              which.max(h)) # highest leverage
names(omit.ind) <- c("infl", "lev")
omit.ind






yobs <- lforest[,"NeedleArea"] #observed values
Mf <- lm(NeedleArea ~ Height + TrunkSize, data = lforest) # all observations
yfitf <- predict(Mf) # fitted values

par(mfrow = c(2,2))
par(mar = c(4,4,0,0)+.1)
cex <- .8
clrs2 <- c("red", "blue")




strikeM2 = data.frame(strike$Strike, strike$Centr, strike$Infl)
pairs(strikeM2, pch = 21, bg = clrs)

strikeMfwd = data.frame(strike$Strike, strike$Country, strike$Infl, strike$Unemp)
pairs(strikeMfwd, pch = 21, bg = clrs)

press1 <- resid(Mfwd)/(1-hatvalues(Mfwd)) # Mfwd
press2 <-  resid(M2)/(1-hatvalues(M2)) # M2
c(PRESS1 = sum(press1^2), PRESS2 = sum(press2^2)) # favors M2