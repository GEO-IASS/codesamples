# Error in y
plot(1, type = 'n', xlim = c(min(sigma_list), max(sigma_list)), ylim = c(0, max(errfit_m)), xlab = "Noise Standard Dev", ylab = "Relative Error", main = "Synthetic Data, Response Error")

cl <- rainbow(CASES)

for (k in 1:CASES) {
  lines(x = sigma_list, y = errfit_m[, k], col = cl[k])
}
legend("topleft", legend = c("OLS", "LASSO", "BLR", "Spike"), col = cl, lwd = 1, cex = 0.5)

# Error in beta
plot(1, type = 'n', xlim = c(min(sigma_list), max(sigma_list)), ylim = c(0, max(errbeta_m)), xlab = "Noise Standard Dev", ylab = "Relative Error", main = "Synthetic Data, Coefficient Error")

cl <- rainbow(CASES)

for (k in 1:CASES) {
  lines(x = sigma_list, y = errbeta_m[, k], col = cl[k])
}
legend("topleft", legend = c("OLS", "LASSO", "BLR", "Spike"), col = cl, lwd = 1, cex = 0.5)

# Error in sparsity: beta
plot(1, type = 'n', xlim = c(min(sigma_list), max(sigma_list)), ylim = c(0, max(sperr_m)), xlab = "Noise Standard Dev", ylab = "Relative Error", main = "Synthetic Data, Sparsity of Coefficient Error")

cl <- rainbow(CASES)

for (k in 1:CASES) {
  lines(x = sigma_list, y = sperr_m[, k], col = cl[k])
}
legend("topleft", legend = c("OLS", "LASSO", "BLR", "Spike"), col = cl, lwd = 1, cex = 0.5)

# Sparsity beta
plot(1, type = 'n', xlim = c(min(sigma_list), max(sigma_list)), ylim = c(0, max(spfit_m)), xlab = "Noise Standard Dev", ylab = "Relative Sparsity", main = "Synthetic Data, Sparsity of Coefficient Vector")

cl <- rainbow(CASES)

for (k in 1:CASES) {
  lines(x = sigma_list, y = spfit_m[, k], col = cl[k])
}
legend("topleft", legend = c("OLS", "LASSO", "BLR", "Spike"), col = cl, lwd = 1, cex = 0.5)

# Posterior Contraction
plot(1, type = 'n', xlim = c(min(sigma_list), max(sigma_list)), ylim = c(0, max(postcontr_m)), xlab = "Noise Standard Dev", ylab = "Contraction", main = "Synthetic Data, Posterior Contraction")

cl <- rainbow(CASES)

for (k in 1:CASES) {
  lines(x = sigma_list, y = postcontr_m[, k], col = cl[k])
}
legend("topleft", legend = c("BLR Mean", "BLR Truth", "Spike Mean", "Spike Truth"), col = cl, lwd = 1, cex = 0.5)
