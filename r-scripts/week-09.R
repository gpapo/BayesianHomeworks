library(tidyverse)
library(compiler)
library(coda)

# Load data
dataset <- read.table("data/msparrownest.dat", col.names = c("nest", "wspan"))

X <- model.matrix(nest ~ wspan, data = dataset)
y <- dataset$nest

n <- NROW(X)
p <- NCOL(X)

# Hyperparameter setting
# pmn_beta <- c(1.1, 0)
pmn_beta <- c(-0.245, 0)
psd_beta <- c(0.3675, 0.049)

invlogit <- function(x) exp(x) / (1 + exp(x))

MHsim <- function(tuning, nsim = 1E4, seed = 2018, verbose = FALSE) {
   varprop <- tuning

   beta <- rep(0, p)
   BETA <- matrix(NA, nsim, p, dimnames = list(NULL, colnames(X)))

   # acceptance counter
   acc <- 0

   set.seed(seed)

   for (iter in 1:nsim) {
      betaprop <- MASS::mvrnorm(1, beta, varprop)
      lhr <- sum(dbinom(y, 1, invlogit(X %*% betaprop), log = TRUE)) +
         sum(dnorm(betaprop, pmn_beta, psd_beta, log = TRUE)) -
         sum(dbinom(y, 1, invlogit(X %*% beta), log = TRUE)) -
         sum(dnorm(beta, pmn_beta, psd_beta, log = TRUE))

      if (log(runif(1)) < lhr) {
         beta <- betaprop
         acc  <- acc + 1
      }

      BETA[iter, ] <- beta
   }

   ESS <- effectiveSize(BETA)
   if (verbose) {
      cat("# nsimulations = ", nsim, "\n# prop. variance:\n", sep = "")
      print(varprop, digits = 2)
      cat("\n# acceptance rate: ", acc / nsim,
          "\n# effective sample size: \n", sep = "")
      print(ESS, digits = 4)
   }
   # [result]
   list(beta = BETA, accept.rate = acc / nsim, eff_samplesize = ESS)
}

nsim <- 1000
varprop <- var(log(y + 1 / 2)) * solve(crossprod(X))
simulation1 <- MHsim(varprop, nsim = nsim, verbose = TRUE)

nsim <- 5 * nsim
varprop <- varprop * 3
simulation2 <- MHsim(varprop, nsim = nsim, verbose = TRUE)

varprop <- varprop * 3
simulation3 <- MHsim(varprop, nsim = nsim, verbose = TRUE)

# Thinning
skeep <- seq(5, nsim, 10)

pdf("week-09_iterations.pdf")
opar <- par(mfrow = c(2, 1))
plot(skeep, simulation3$beta[skeep, 1],
     type = "l", xlab = "Iteration", ylab = expression(alpha))
plot(skeep, simulation3$beta[skeep, 2],
     type = "l", xlab = "Iteration", ylab = expression(beta))
par(opar)
graphics.off()

# Comparison
alpha <- seq(-10, 10, by = 0.1)
beta  <- seq(-2, 2, by = 0.1)

pdf("week-09_priors_posteriors.pdf", height = 6, width = 10)
opar <- par(mfrow = c(1, 2))

plot(alpha, dnorm(alpha, pmn_beta[1], psd_beta[1]),
     lwd = 1, lty = 2, ylim = c(0, 0.25),
     type = "l", xlab = expression(alpha), ylab = "density")
lines(density(simulation3$beta[, 1]), lty = 1)
legend("topright", legend = c("prior", "posterior"), cex = 1,
       lty = c(2, 1), seg.len = 1.5, bty = "n")

plot(beta, dnorm(beta, pmn_beta[2], psd_beta[2]),
     lwd = 1, lty = 2, ylim = c(0, 3.2),
     type = "l", xlab = expression(beta), ylab = "density")
lines(density(simulation3$beta[, 2], lty = 1))
legend("topright", legend = c("prior", "posterior"), cex = 1,
       lty = c(2, 1), seg.len = 1.5, bty = "n")

par(opar)
graphics.off()

ysim <- invlogit(tcrossprod(simulation3$beta, X))
CI <- apply(ysim, 2, quantile, prob = c(0.025, 0.5, 0.9975))

pdf("week-09_predictions.pdf", height = 4)
plot(c(10, 15), range(c(0, CI)), type = "n", xlab = "wingspan", ylab = "f")
lines(CI[1, ], col = "blue", lwd = 2)
lines(CI[2, ], col = "red",  lty = 2)
lines(CI[3, ], col = "blue", lwd = 2)
graphics.off()
