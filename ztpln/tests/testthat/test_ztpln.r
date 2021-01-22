library(ztpln)

seed <- 123

# K = 1 --------------------------
set.seed(seed)
K <- 1
mu <- runif(K, 0, 15)
sig <- runif(K, 0, 2)
n <- 100

mu
sig

y <- rztpln(n = 100, mu, sig)
y

lik1 <- sum(dztpln(y, mu, sig, log = TRUE))
#lik2 <- sum(log(poilog::dpoilog(y, mu, sig) / (1 - poilog::dpoilog(0, mu, sig))))

k <- 5
mu <- 3
sig <- 2
y <- rztpln(n = 100, mu, sig)
y

lik1 <- sum(dztpln(k, mu, sig, log = TRUE))
#lik2 <- sum(log(poilog::dpoilog(k, mu, sig) / (1 - poilog::dpoilog(0, mu, sig))))

#poilog::dpoilog(0, mu, sig)
#do_dpln(0, mu, sig^2)

#lik1 - lik2

ztplnMLE(y)


# K = 2 --------------------------
set.seed(seed)
K <- 2
mu <- runif(K, 0, 15)
sig <- runif(K, 0, 2)
n <- 100

# random draw from dirichlet  distribution
theta0 <- rgamma(K, sample(1:3, K, replace = TRUE))
theta <- theta0 / sum(theta0)

mu
sig
theta

y <- rztplnm(n = 100, mu, sig, theta)
y

lik1 <- sum(dztplnm(y, mu, sig, theta, log = TRUE))

tmp <- NULL
for (k in 1:K) tmp <- cbind(tmp, theta[k] * dztpln(y, mu[k], sig[k]))
lik2 <- sum(log(apply(tmp, 1, sum)))

lik1 - lik2

ztplnmMLE(y, K = 2)
ztplnmMLE(y, K = 2, type1 = FALSE)


# K = 3 --------------------------
set.seed(seed)
K <- 2
mu <- runif(K, 0, 15)
sig <- runif(K, 0, 2)
n <- 100

# random draw from dirichlet  distribution
theta0 <- rgamma(K, sample(1:3, K, replace = TRUE))
theta <- theta0 / sum(theta0)

mu
sig
theta

y <- rztplnm(n = 100, mu, sig, theta)
y

lik1 <- sum(dztplnm(y, mu, sig, theta, log = TRUE))

tmp <- NULL
for (k in 1:K) tmp <- cbind(tmp, theta[k] * dztpln(y, mu[k], sig[k]))
lik2 <- sum(log(apply(tmp, 1, sum)))

lik1 - lik2

ztplnmMLE(y, K = 3)
ztplnmMLE(y, K = 3, type1 = FALSE)
