#
# Description:
# Exercise 9.2, Hoff's book
#
# Author: Giovanni Papini
# Date:   2018-11-21
#

library(tidyverse)
library(compiler)

# PART A ####
set.seed(1)

# loading data
dataset <- read_csv("data/azdiabetes.dat", col_names = TRUE)
N <- nrow(dataset)

y <- dataset$glu
X <- model.matrix(glu ~ ., data = dataset)

# rescale when numeric: (y - mean) / sd
y <- scale(y)
X <- apply(X, 2, function(x) if (length(unique(x)) > 2) scale(x) else x)

# OLS modeling
lsfit <- lm(y ~ X - 1)

# G prior
lmgprior <- function(y, X, g = NROW(X),
                nu_0 = 1, s2_0 = try(summary(lm(y ~ X - 1))$sigma ^ 2, silent = TRUE),
                S = 5000) {
   n <- NROW(X)
   p <- NCOL(X)

   H_g   <- (g / (g + 1)) * X %*% solve(t(X) %*% X) %*% t(X)
   SSR_g <- c(t(y) %*% (diag(n) - H_g) %*% y)
   s2    <- 1 / rgamma(S, (nu_0 + n) / 2, (nu_0 * s2_0 + SSR_g) / 2)

   V_b <- g * solve(t(X) %*% X) / (g + 1)
   E_b <- V_b %*% t(X) %*% y # Expected value of the posterior

   # Use the s2 sampled from the inv-gamma
   E <- matrix(rnorm(S * p, 0, sqrt(s2)), S, p)

   # Sampling from a normal centered in E_b and with variance s2 * V_b
   # beta = E_b + E * chol(V_b), E ~ N(0, s2 * I)
   beta <- t(t(E %*% chol(V_b)) + c(E_b))

   # Return
   list(beta = beta, s2 = s2)
}

lmgpriorfit <- lmgprior(y, X, g = N, nu_0 = 2, s2_0 = 1)

beta_gprior <- apply(lmgpriorfit$beta, 2,
                     function(x) quantile(x, c(0.025, 0.975)))
#xtable::xtable(beta_gprior)

# PART B ####
set.seed(1)

backselect_tcrit <- function(y, X, tcrit = qnorm(0.975)) {
   Xc <- X
   remaining <- 1:NCOL(X)

   repeat {
      if (!is_empty(remaining)) {
         fit   <- lm(y ~ Xc - 1)
         summ  <- summary(fit)
         tstat <- summ$coef[, 3]
         s2    <- summ$sigma ^ 2
      } else break
      jpmax <- which.min(abs(tstat))
      Xc <- Xc[, -jpmax, drop = FALSE]
      remaining <- remaining[-jpmax]

      if ( (all(abs(tstat) > tcrit) || is_empty(remaining)) ) break
   }

   # Return
   list(remaining = remaining, removed = setdiff(1:NCOL(X), remaining))
}

loglikpy_X <- function(y, X, g = NCOL(X),
                       nu_0 = 1,
                       s2_0 = try(summary(lm(y ~ X - 1))$sigma ^ 2,
                                  silent = TRUE)) {
   n <- NROW(X)
   p <- NCOL(X)

   if (p == 0) s2_0 <- mean(y ^ 2)
   H_0 <- 0
   if (p > 0) H_0 <- g / (g + 1) * X %*% solve(t(X) %*% X) %*% t(X)
   SSR_g <- c(t(y) %*% (Id(n) - H_0) %*% y)

   # Return
   - 0.5 * n * log(2 * pi) + lgamma(0.5 * (nu_0 + n)) - lgamma(0.5 + nu_0) +
      - 0.5 * p * log(1 + g) + 0.5 * nu_0 * log(0.5 * nu_0 * s2_0) +
      - 0.5 * (nu_0 + n) * log(0.5 * (nu_0 * s2_0 + SSR_g))
}

# model selection (backward elimination)
variables <- backselect_tcrit(y, X, tcrit = 1.65)
Xred <- X[, variables$remaining]
reducedlmfit <- lm(y ~ Xred - 1)

# model averaging
p <- NCOL(X)
S <- 500
B <- matrix(NA, S, p, dimnames = list(NULL, colnames(X)))
Z <- matrix(NA, S, p, dimnames = list(NULL, colnames(X)))
z <- rep(TRUE, p)

lpy_0 <- loglikpy_X(y, X[, z])

for (iter in 1:S) {
   for (j in sample(1:p)) {
      zp    <- z # include/exclude a sampled regressor
      zp[j] <- !zp[j]
      lpy_p <- loglikpy_X(y, X[, zp])
      r     <- (lpy_p - lpy_0) * (-1) ^ (!zp[j])

      # Note: Prior on z gives same prior probability to all models
      z[j] <- rbinom(1, 1, 1 / (1 + exp(-r))) == 1
      if (z[j] == zp[j]) lpy_0 <- lpy_p
   }

   beta <- z
   if (sum(z) > 0) {
      beta[z] <- lmgprior(y, X[, z, drop = FALSE], S = 1)$beta
   }

   Z[iter, ] <- z
   B[iter, ] <- beta
}

beta_bma <- apply(B, 2, quantile, c(0.025, 0.975), na.rm = TRUE)
prob_include <- apply(Z, 2, mean)

#xtable::xtable(beta_bma)
#xtable::xtable(t(prob_include))
