# Florica Constantine 
# High Dimensional Regression Exploration 
# Project for STAT 632

library(glmnet) # LASSO
library(MCMCpack) # Bayesian Linear Regression
library(BoomSpikeSlab) # Bayesian Linear Regression with Spike and Slab Prior
library(hdi) # High Dimensional Inference: For Data

# Real Data
data(riboflavin)
attach(riboflavin)
X <- data.matrix(x)
y <- riboflavin$y

# Make dataframe
yxdf <- as.data.frame(cbind(y, X))
names(yxdf)[1] <- "y"

# Parameters
p = dim(X)[2] # Variables
n = length(y) # Samples

thresh = log(p) / p # For determining sparsity

# Number of Cases
CASES <- 4

# Outputs

# Error in Fit Relative to y
errfit <- vector(mode = "numeric", length = CASES)
# Sparsity of betahat
spfit <- vector(mode = "numeric", length = CASES)
# Posterior Contraction of betahat Samples
postcontr <- vector(mode = "numeric", length = CASES / 2)

# Least Squares

# Coefficients
betahat <- ginv(X) %*% y
# Error in Response
errfit[1] = sqrt(sum((y - (X %*% betahat))^2))
# Sparsity of Coefficients 
spfit[1] = length(which(abs(betahat) > thresh))

# LASSO

# Find LASSO object
cvfit <- cv.glmnet(x = X, y = y)
# Coefficients
betahat <- coef(cvfit, s = "lambda.1se") # ``Best''
betahat <- matrix(betahat)[-1]
# Find Error in Fit
errfit[2] <- sqrt(sum((y - (X %*% betahat))^2))
# Sparsity of Coefficients 
spfit[2] = length(which(abs(betahat) > thresh))

# Bayesian Regression

# Fit Bayesian Regression Model
blrfull <- MCMCSVDreg(y ~ ., data = yxdf, burnin = 500, mcmc = 2000, thin = 1, verbose = 0, tau2.start = 1.0, g0 = 0.0, a0 = 0.001, b0 = 0.001, c0 = 2.0, d0 = 2.0, w0 = 0.0, beta.samp = TRUE, intercept = FALSE)
blr <- blrfull[,(dim(blrfull)[2] - p + 1):dim(blrfull)[2]]
# Use Posterior Mean as Fit
betahat <- apply(blr, 2, mean)
# Find Error in Fit
errfit[3] <- sqrt(sum((y - (X %*% betahat))^2))
# Sparsity of Coefficients
spfit[3] = length(which(abs(betahat) > thresh))
# Posterior Contraction Around Mean
postcontr[1] <- mean(sqrt(apply(apply(blr, 1, '-', betahat)^2, 1, sum)))

# Bayesian Spike-and-Slab Regression

# Fit Model
spl <- lm.spike(y ~ ., niter = 2500, data = yxdf, ping = 10, bma.method = "ODA", error.distribution = "gaussian")
# Coefficients: Posterior Mean
bt <- spl$beta[501:2500, 2:(p + 1)]
betahat <- apply(bt, 2, mean)
# Find Error in Fit
errfit[4] <- sqrt(sum((y- (X %*% betahat))^2))
# Sparsity of Coefficients
spfit[4] = length(which(abs(betahat) > thresh))
# Posterior Contraction Around Mean
postcontr[2] <- mean(sqrt(apply(apply(bt, 1, '-', betahat)^2, 1, sum)))

errfit = errfit / sqrt(sum(y^2))
