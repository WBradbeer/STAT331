
# load file
strike <- read.csv("strikes_clean.csv")
# encode categoricals as Factors
strike$Centr <- as.factor(strike$Centr)
strike$Country <- as.factor(strike$Country)

# Investigate to look for patterns
summary(strike)
pairs(strike)

# remove zero values
strike <- strike[strike$Strike != 0,]


# Reason that these should affect strike
M1 = lm(Strike ~ Infl*Centr*Unemp, data = strike)
# anova shows unemp doesn't add much 
anova(M1)
# remove Unemp from M1 to get current model that works well 
M2 = lm(Strike ~ Infl*Centr, data = strike)
# multiplicative model
M3 = lm(log(Strike) ~ Infl*Centr, data = strike)

M0.log = lm(log(Strike) ~ 1, data = strike)

Mfull.log = lm(log(Strike) ~ (. - Country)^2, data = strike)
# basic model
M0 <- lm(Strike ~ 1, data = strike)
# full model without country - interaction effects
Mfull <- lm(Strike ~ (. - Country - Year)^2 , data = strike)
# forward step model from intercept to full
Mfwd <- step(M0, scope=list(lower=M0, upper=Mfull), direction='forward', trace=FALSE)
# backward step model from full to intercept
Mback <- step(Mfull, scope=list(lower=M0, upper=Mfull), direction='backward', trace=FALSE)
# step model starting from model grenerated from exploration
Mstep <- step(M2, scope=list(lower=M0, upper=Mfull), direction='both', trace=FALSE)

Mfwd.log <- step(M0.log, scope=list(lower=M0.log, upper=Mfull.log), direction='forward', trace=FALSE)
# backward step model from full to intercept
Mback.log <- step(Mfull.log, scope=list(lower=M0.log, upper=Mfull.log), direction='backward', trace=FALSE)
# step model starting from model grenerated from exploration
Mstep.log <- step(M3, scope=list(lower=M0.log, upper=Mfull.log), direction='both', trace=FALSE)


models <- list(M0, M1, M2, Mfwd, Mback, Mstep, M3, Mfwd.log, Mback.log, Mstep.log)


n <- length(strike$Strike)
samp.size <- 400
results <- rep(0, length(models))

for(i in 1:2000) {
  # seperate training and test data
  samp <- sample(n, samp.size)
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

# model diagnostics
res <- resid(Mback)
H <-  model.matrix(Mback)
H <- H %*% solve(crossprod(H), t(H))
h <- diag(H)
# studentized residuals
res.stu <- res/sqrt(1-h)

# standardized student residuals
res.stand <- res.stu/sigma(Mback)

# diagnostic plots
plot(predict(Mback), res.stu) # seems like multiplicative error
hist(res.stu,breaks=100)

# try multiplicative model
strike.nonzero <- strike[strike$Strike != 0,]
M3 = lm(log(Strike) ~ Infl*Centr, data = strike.nonzero)

