
# load file
strike <- read.csv("strikes_clean.csv")
# encode categoricals as Factors
strike$Centr <- as.factor(strike$Centr)
strike$Country <- as.factor(strike$Country)

# Investigate to look for patterns
summary(strike)
pairs(strike)

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
