
library(tidyverse)
library(mvnfast)
library(compiler)

invlogit <- function(x) 1 / (1 + exp(-x))

# Load data
mathstandard <- read.table("data/mathstandard2.dat", header = TRUE)

# b)
GLMfits <- by(mathstandard, with(mathstandard, county), function(data)
   glm(metstandard ~ percentms, family = binomial(), data = data))

Beta     <- t(sapply(GLMfits, coefficients))
n        <- with(mathstandard, table(county))
topn     <- n %>% keep(~ . >= 10)
ML_theta <- Beta[names(topn), ] %>% colMeans
ML_Sigma <- Beta[names(topn), ] %>% cov

palette(grey.colors(max(n), 0.8, 0))

pdf("week-10_groupedglm.pdf", height = 5)

plot(c(0, 100), c(0, 1), type = 'n', bty = 'n',
     xlab = "teachers with MS (%)", ylab = "P(Y = 1)")
for (i in 1:NROW(Beta)) {
   curve(invlogit(Beta[i, 1] + Beta[i, 2] * x), add = TRUE, col = n[i])
}

curve(invlogit(ML_theta[1] + ML_theta[2] * x), add = TRUE,
      col = "red", lty = 2, lwd = 4)

graphics.off()

pdf("week-10_frequencies.pdf", height = 4)
tibble(county = names(n), freq = n, min5 = factor(n < 5)) %>%
   ggplot(aes(county, freq, fill = min5)) +
   geom_col(show.legend = FALSE) +
   theme_minimal() + scale_fill_manual(values = c("#1E90FF", "#FF0000")) +
   theme(axis.text.x = element_text(angle = 90)) +
   ylab("Frequency") + xlab("County")
graphics.off()

# c)
Y <- with(mathstandard, tapply(metstandard, county, identity))
X <- with(mathstandard, split(mathstandard, county)) %>%
   map(~ model.matrix(metstandard ~ percentms, data = .))

p <- 2
m <- length(X)

simulateMH <- function(niter, init, seed = NULL) {
   set.seed(seed)

   # Hyperparameters
   mu0   <- init$mu
   S0    <- init$S
   eta0  <- init$eta
   Beta  <- init$Beta
   invL0 <- invSigma <- solve(S0)

   # Init storing variables
   THETA <- matrix(NA, niter, p)
   BETA  <- array(NA, c(m, p, niter))
   SIGMA <- array(NA, c(p, p, niter))

   colnames(THETA) <- names(mu0)
   dimnames(BETA)  <- dimnames(Beta)
   dimnames(SIGMA) <- dimnames(S0)

   for (itr in 1:niter) {
      # update theta
      Lm    <- solve(invL0 + m * invSigma)
      mu_m  <- Lm %*% (invL0 %*% mu0 + invSigma %*% colSums(Beta))
      theta <- rmvn(1, mu_m, Lm) %>% drop()

      # update Sigma
      tmp      <- sweep(Beta, 2, theta)
      invSigma <- rWishart(1, eta0 + m, solve(S0 + crossprod(tmp))) %>% drop()
      Sigma    <- solve(invSigma)

      # update Beta
      Beta0      <- tapply(Beta, row(Beta), t)
      Beta1      <- tapply(Beta, row(Beta), rmvn, n = 1, sigma = Sigma / 2)
      lpBeta0    <- sapply(Beta0, dmvn, theta, Sigma, log = TRUE)
      lpBeta1    <- sapply(Beta1, dmvn, theta, Sigma, log = TRUE)
      pi.Beta0   <- mapply(tcrossprod, X, Beta0, SIMPLIFY = FALSE) %>% lapply(invlogit)
      pi.Beta1   <- mapply(tcrossprod, X, Beta1, SIMPLIFY = FALSE) %>% lapply(invlogit)
      lpyi.Beta0 <- mapply(dbinom, Y, pi.Beta0, MoreArgs = list(size = 1, log = TRUE),
                           SIMPLIFY = FALSE)
      lpyi.Beta1 <- mapply(dbinom, Y, pi.Beta1, MoreArgs = list(size = 1, log = TRUE),
                           SIMPLIFY = FALSE)

      lr <- sapply(lpyi.Beta1, sum) - sapply(lpyi.Beta0, sum) + lpBeta1 - lpBeta0
      lu <- log(runif(m))
      Beta[lu < lr, ] <- do.call(rbind, Beta1[lu < lr])

      THETA[itr, ]   <- theta
      BETA[, , itr]  <- Beta
      SIGMA[, , itr] <- Sigma
   }

   # [return]
   list(THETA = THETA, SIGMA = SIGMA, BETA = BETA)
}

init  <- list(mu = ML_theta, S = ML_Sigma, eta = 4, Beta = replace_na(Beta, 0))
niter <- 10000

mcmcsim <- simulateMH(niter, init)

# d)

B <- 500 # burn-in

pdf("week-10_mcmcplot.pdf", height = 9)

mcmcsim %$%
   tibble(`theta[intercept]`   = THETA[, 1],
          `theta[percentms]`   = THETA[, 2],
          `sigma[intercept]^2` = SIGMA[1, 1, ],
          `sigma[int * pms]`   = SIGMA[1, 2, ],
          `sigma[percentms]^2` = SIGMA[2, 2, ],
          iteration            = seq_len(niter),
          status               = factor(iteration < B, c(TRUE, FALSE),
                                        c("`Burn-in`", "Warmed"))) %>%
   gather(param, value, -iteration, -status) %>%
   ggplot(aes(iteration, value)) + geom_line() +
   facet_wrap(~ param + status, nc = 2, scale = "free",
              labeller = label_parsed) +
   theme_minimal() + ylab("") + xlab("")

graphics.off()

E_Beta  <- mcmcsim %$% apply(BETA[, , B:niter], 1:2, mean)
E_Sigma <- mcmcsim %$% apply(SIGMA[, , B:niter], 1:2, mean)
E_theta <- mcmcsim %$% colMeans(THETA[B:niter, ])
V_theta <- mcmcsim %$% cov(THETA[B:niter, ])

pdf("week-10_groupedposterior.pdf", height = 5)

plot(c(0, 100), c(0, 1), type = 'n', bty = 'n',
     xlab = "teachers with MS (%)", ylab = "P(Y = 1)")
for (i in 1:NROW(E_Beta)) {
   curve(invlogit(E_Beta[i, 1] + E_Beta[i, 2] * x), add = TRUE, col = n[i])
}

curve(1 / (1 + exp(-E_theta[1] - E_theta[2] * x)), add = TRUE,
      col = "red", lty = 2, lwd = 3)
curve(1 / (1 + exp(-ML_theta[1] - ML_theta[2] * x)), add = TRUE,
      col = "brown", lty = 3, lwd = 3)

graphics.off()

# f)

pdf("week-10_prior_posterior-1.pdf", height = 5, width = 8)
par(mfrow = c(1, 2))

curve(dnorm(x, ML_theta[1], sqrt(ML_Sigma[1, 1])), -10, 5,
      ylim = c(0, 0.8), bty = 'n',
      main = "", ylab = "", xlab = expression(theta[intercept]))
mcmcsim %$% density(THETA[B:niter, 1]) %>% lines(col = c("#D60000"))

legend("topleft", legend = c("prior", "posterior"),
       col = c("black", "#D60000"), lty = 1, bty = 'n')

curve(dnorm(x, ML_theta[2], sqrt(ML_Sigma[2, 2])), -0.1, 0.15,
      ylim = c(0, 50), bty = 'n',
      main = "", ylab = "", xlab = expression(theta[percentms]))
mcmcsim %$% density(THETA[B:niter, 2]) %>% lines(col = c("#D60000"))

graphics.off()
