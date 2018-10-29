#
# Description:
# Week 06 homeworks for Bayesian Inference course
#
# Author: Giovanni Papini
# Date:   2018-10-28
#

rm(list = ls())

library(tidyverse)

# Load data
yobs <- list(
 A = scan("data/menchild30bach.dat"),
 B = scan("data/menchild30nobach.dat")
)

# Hyperparameters
a_theta <- 2
b_theta <- 1
gammapriors <- 2 ^ (3:7)

# Sample statistics
ytot <- map(yobs, sum)
n <- map(yobs, length)

# Simulation main function
library(compiler) # to accelerate
gibbs_sim <- function(a_gamma, b_gamma = a_gamma, theta0 = 1,
                      gamma0 = 1, nsim = 5000, seed = 10) {
   set.seed(seed)
   # Initialize
   result <- matrix(NA, nsim, 2, dimnames = list(NULL, c("theta", "gamma")))
   result[1, ] <- c(theta0, gamma0)
   # Main loop
   for (r in 2:nsim) {
      result[r, "theta"] <- rgamma(1, a_theta + ytot$A + ytot$B,
                                   b_theta + n$A + n$B * result[r - 1, "gamma"])
      result[r, "gamma"] <- rgamma(1, a_gamma + ytot$B,
                                   b_gamma + n$B * result[r, "theta"])
   }
   # Return
   return(as_tibble(result))
}

# Simulation
simulations <- map(gammapriors, gibbs_sim, nsim = 1E4)
statistics <- map_dbl(simulations, ~ with(., mean(theta * gamma - theta)))
# 0.3720311 0.3344181 0.2713526 0.2000735 0.1327382

# Plots
pdf("week-06_hist.pdf")
opar <- par(mfrow = n2mfrow(6))

for (j in seq_along(gammapriors)) {
   theta <- pluck(simulations, j, "theta")
   gamma <- pluck(simulations, j, "gamma")
   hist(theta, prob = TRUE, col = "lightblue", ylim = c(0, 5), xlim = c(0, 2.5),
        ylab = "posterior density", xlab = "mean number of children",
        main = substitute(paste(a[gamma], " = ", b[gamma], " = " , prior),
                          list(prior = gammapriors[j])))
   lines(density(theta), col = "blue")
   hist(gamma * theta, prob = TRUE, col = "orange", add = TRUE)
   lines(density(gamma * theta), col = "red")
   legend("topleft", c(expression(theta[A]), expression(theta[B])),
          col = c("blue", "red"), lty = 1, cex = 0.7)
}

par(opar)
graphics.off()

pdf("week-06_prior.pdf")

plot(NULL, xlim = c(0, 3), ylim = c(0, 4.5),
     xlab = expression(gamma), ylab = expression(p(gamma)),
     main = "Gamma prior")
for (i in seq_along(gammapriors)) {
   curve(dgamma(x, gammapriors[i], gammapriors[i]), add = TRUE, col = i)
}
legend_txt <- parse(text = paste("a[gamma] ==", gammapriors))
legend(x = 1.5, y = 3, legend_txt, col = seq_along(gammapriors), lty = 1)
mtext(expression(b[gamma] == a[gamma]))

graphics.off()

# Diagnostics
library(coda)
map(simulations, effectiveSize)

# [[1]]
# theta    gamma
# 1175.254 1198.183
#
# [[2]]
# theta    gamma
# 1471.870 1353.061
#
# [[3]]
# theta    gamma
# 1840.738 1695.119
#
# [[4]]
# theta    gamma
# 2424.931 2238.538
#
# [[5]]
# theta    gamma
# 3259.088 2942.661
