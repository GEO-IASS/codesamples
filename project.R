# Florica Constantine 
# High Dimensional Regression Exploration 
# Project for STAT 632

library(glmnet) # LASSO
library(MCMCpack) # Bayesian Linear Regression
library(BoomSpikeSlab) # Bayesian Linear Regression with Spike and Slab Prior

# Synthetic Data 

# Parameters 
n <- 50 # Samples
p <- 100 # Dimension
sp <- round(0.05 * p) # Sparsity
thresh <- log(p) / p # To assess sparsity

# Design Matrix
X <- matrix(rnorm(p * n, mean = 0.0, sd = 1.0), nrow = n, ncol = p)
X[, p] <- 1.0 # Last Column is Intercept

# Coefficient Vector 
beta <- rep(0, p)
beta[1:sp] <- 1.0 
delta <- which(abs(beta) > thresh) # Sparsity pattern 

# Response without Noise
y_clean <- X %*% beta

# Noise Levels
sigma_list <- seq(0, 5.0 * sqrt(sum(beta^2)), length = 10)
s_l <- length(sigma_list)

# Trials
trials <- 20

# Number of Cases
CASES <- 4

# Outputs

# Error in Fit Relative to y_clean
errfit <- array(0, c(s_l, trials, CASES))
# Error in Fit Relative to beta
errbeta <- array(0, c(s_l, trials, CASES))
# Error in Sparsity Pattern of betahat
sperr <- array(0, c(s_l, trials, CASES))
# Sparsity of betahat
spfit <- array(0, c(s_l, trials, CASES))
# Posterior Contraction of betahat Samples
postcontr <- array(0, c(s_l, trials, CASES))

# Simulation

for (ss in 1:s_l) { # Noise Level
	sigma_e <- sigma_list[ss]

	for (tr in 1:trials) { # Trials
	  cat(tr, ss)
	  
		# Noisy Data
		y <- y_clean + c(rnorm(n, mean = 0.0, sd = sigma_e))
	
		# Make dataframe
		yxdf <- as.data.frame(cbind(y, X))
		names(yxdf)[1] <- "y"
		
		# Least Squares
		
		# Coefficients
		betahat <- ginv(X) %*% y;
		# Find Error in Fit
		errfit[ss, tr, 1] <- sqrt(sum((y_clean - (X %*% betahat))^2))
		# Find Error in beta
		errbeta[ss, tr, 1] <- sqrt(sum((betahat - beta)^2))
		# Find Sparsity of Fit
		delta_hat <- which(betahat > thresh)
		# Error in Sparsity Pattern
		sperr[ss, tr, 1] <- length(union(delta, delta_hat)) - sp
		# Sparsity of Fit
		spfit[ss, tr, 1] <- length(delta_hat)

		# LASSO

		# Find LASSO object
		cvfit <- cv.glmnet(x = X, y = y)
		# Coefficients
		betahat <- coef(cvfit, s = "lambda.1se") # ``Best''
		betahat <- matrix(betahat)[-1]
		# Find Error in Fit
		errfit[ss, tr, 2] <- sqrt(sum((y_clean - (X %*% betahat))^2))
		# Find Error in beta
		errbeta[ss, tr, 2] <- sqrt(sum((betahat - beta)^2))
		# Find Sparsity of Fit
		delta_hat <- which(abs(betahat) > thresh)
		# Error in Sparsity Pattern
		sperr[ss, tr, 2] <- length(union(delta, delta_hat)) - sp
		# Sparsity of Fit
		spfit[ss, tr, 2] <- length(delta_hat)

		# Bayesian Regression
		
		# Fit Bayesian Regression Model
		blrfull <- MCMCSVDreg(y ~ ., data = yxdf, burnin = 500, mcmc = 2000, thin = 1, verbose = 250, tau2.start = 1.0, g0 = 0.0, a0 = 0.001, b0 = 0.001, c0 = 2.0, d0 = 2.0, w0 = 0.0, beta.samp = TRUE, intercept = FALSE)
		blr <- blrfull[,(dim(blrfull)[2] - p + 1):dim(blrfull)[2]]
		# Use Posterior Mean as Fit
		betahat <- apply(blr, 2, mean)
		# Find Error in Fit
		errfit[ss, tr, 3] <- sqrt(sum((y_clean - (X %*% betahat))^2))
		# Find Error in beta
		errbeta[ss, tr, 3] <- sqrt(sum((betahat - beta)^2))
		# Find Sparsity of Fit
		delta_hat <- which(abs(betahat) > thresh)
		# Error in Sparsity Pattern
		sperr[ss, tr, 3] <- length(union(delta, delta_hat)) - sp
		# Sparsity of Fit
		spfit[ss, tr, 3] <- length(delta_hat)
		# Posterior Contraction Around Mean
		postcontr[ss, tr, 1] <- mean(sqrt(apply(apply(blr, 1, '-', betahat)^2, 1, sum)))
		# Posterior Contraction Around Truth
		postcontr[ss, tr, 2] <- mean(sqrt(apply(apply(blr, 1, '-', beta)^2, 1, sum)))

		# Bayesian Spike-and-Slab Regression

		# Fit Model
		spl <- lm.spike(y ~ ., niter = 2500, data = yxdf, ping = 250, bma.method = "ODA", error.distribution = "gaussian")
		# Coefficients: Posterior Mean
		bt <- spl$beta[501:2500, 2:(p + 1)]
		betahat <- apply(bt, 2, mean)
		# Find Error in Fit
		errfit[ss, tr, 4] <- sqrt(sum((y_clean - (X %*% betahat))^2))
		# Find Error in beta
		errbeta[ss, tr, 4] <- sqrt(sum((betahat - beta)^2))
		# Find Sparsity of Fit
		delta_hat <- which(abs(betahat) > thresh)
		# Error in Sparsity Pattern
		sperr[ss, tr, 4] <- length(union(delta, delta_hat)) - sp
		# Sparsity of Fit
		spfit[ss, tr, 4] <- length(delta_hat)
		# Posterior Contraction Around Mean
		postcontr[ss, tr, 3] <- mean(sqrt(apply(apply(bt, 1, '-', betahat)^2, 1, sum)))
		# Posterior Contraction Around Truth
		postcontr[ss, tr, 4] <- mean(sqrt(apply(apply(bt, 1, '-', beta)^2, 1, sum)))
	}
}

# Relativize Errors
errfit_m <- apply(errfit / sqrt(sum(y_clean^2)), c(1, 3), mean)
errbeta_m <- apply(errbeta / sqrt(sum(beta^2)), c(1, 3), mean)
sperr_m <- apply(sperr / sp, c(1, 3), mean)
spfit_m <- apply(spfit / sp, c(1, 3), mean)
postcontr_m <- apply(postcontr / sqrt(sum(beta^2)), c(1, 3), mean)
