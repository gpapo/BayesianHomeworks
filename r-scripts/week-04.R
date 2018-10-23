
# Exercise 1 --------------------------------------------------------------

set.seed(10)
library(compiler)
count_switch <- function(x) x %>% diff() %>% abs() %>% sum()

y <- c(1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0)
t_obs <- sum(abs(diff(y)))

# posterior hyperparameters
alpha <- 1 + sum(y)
beta  <- 1 + length(y) - sum(y)

niter <- 10000
ssize <- 20

theta <- rbeta(niter, alpha, beta)
t_rep <- map_dbl(theta, ~ rbinom(ssize, 1, .) %>% count_switch())

# result
pval <- 2 * mean(t_rep <= t_obs)
cat("Bayesian p-val:", pval)

# plots
hist(t_rep)
abline(v = t_obs, col = 2, lty = 2)

# Exercise 2 --------------------------------------------------------------
library(compiler) # for matter of loop performance

envelop <- function(n, f, g, randomg, start, lower = 0, upper = 1 - lower, ...) {
   result <- double(n)
   optimization <- nlminb(objective = function(x) -f(x) / g(x), start, lower = lower, upper = upper, ...)
   M <- -optimization$obj
   ind <- 1
   while (ind <= n) {
      x <- randomg(1)
      u <- runif(1, 0, M * g(x))
      if (u < f(x)) {
         result[ind] <- x
         ind <- ind + 1
      }
   }
   return(result)
}

# custom function
ssize <- 10000
f <- dexp
g <- partial(dbeta, shape1 = 0.5, shape2 = 0.5)
randomg <- partial(rbeta, shape1 = 0.5, shape2 = 0.5)

y <- envelop(ssize, f, g, randomg, start = 0.5)

plot(density(y))
curve(f, add = TRUE, col = 2)

# curve plot
curve(f(x) / g(x), 0, 2, col = "blue", lty = 2, ylim = c(0, 1.5),
      xlab = "x", ylab = "y")
curve(f, col = "red", add = TRUE)
curve(g, col = "green", add = TRUE)
M <- optimise(function(x) -f(x) / g(x), c(0, 1))
abline(h = - M$obj, lty = 3, col = "darkblue")
text(M$min, -M$obj + 0.05, "M")
segments(0.6, 0.0, 0.6, f(0.6), lwd = 2, col = "darkred")
segments(0.6, f(0.6), 0.6, -M$obj * g(0.6), lwd = 2)
segments(-10, -M$obj * g(0.6), 0.6, -M$obj * g(0.6), lty = 3)
