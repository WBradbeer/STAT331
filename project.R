
# load file
strike <- read.csv("strikes_clean.csv")
# encode categoricals as Factors
strike$Centr <- as.factor(strike$Centr)
strike$Country <- as.factor(strike$Country)

# current model that works well 
M = lm(Strike ~ Infl*Centr, data = strike)
