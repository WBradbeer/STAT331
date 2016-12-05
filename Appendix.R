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

# Section 3.1: Residual Plots

# First we find the residual plot of M_a
res_a <- resid(M_a)
H_a <-  model.matrix(M_a)
H_a <- H_a %*% solve(crossprod(H_a), t(H_a))
h_a <- diag(H_a)
# Studentized residuals
res.stu_a <- res_a/sqrt(1-h_a)

# Standardized student residuals
res.stand_a <- res.stu_a/sigma(M_a)

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
res.stand_b <- res.stu_b/sigma(M_b)

cex <- 0.8 
par(mar = c(4,4,0.1,0.1))
plot(predict(M_b), res_b, pch=21, bg="black", cex=cex, cex.axis=cex,  
     xlab="Predicted # of strike-hours", ylab="Residual strike-hours")
points(predict(M_b), res.stu_b, pch=21, bg="red", cex=cex)
legend(x="topleft", c("Residuals", "Studentized Residuals"), pch=21,
       pt.bg = c("black", "red"), pt.cex=cex, cex=cex)



# Next, we will plot the histogram of studentized residuals for M_a
par(mar = c(4,4,0.1,0.1))
hist.plot_a <- hist(res.stand_a, breaks=75, freq = FALSE, cex.axis=cex, main="",
                   xlab="Standardized, studentized residual strike activity")
curve(dnorm(x), col='red', add=TRUE)

# Next, we will plot the histogram of studentized residuals for M_b
par(mar = c(4,4,0.1,0.1))
hist.plot_b <- hist(res.stand_b, breaks=100, freq = FALSE, cex.axis=cex, main = "",
                   xlab="Standardized, studentized residual strike activity")
curve(dnorm(x), col='red', add=TRUE)