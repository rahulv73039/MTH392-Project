rm(list = ls())

#---------------
# CDF, PDF, Quantile and simulation function of Asymmetric Laplace (AL)

pasymlp <- function(x, del.l, del.u){
  coeff <- del.l * as.numeric(x <= 0) + del.u * as.numeric(x > 0)
  out <- coeff * exp(-abs(x) / coeff) / (del.l + del.u)
  out <- ifelse(x <= 0, out, 1 - out)
  return(out)
}

dasymlp <- function(x, del.l, del.u){
  coeff <- del.l * as.numeric(x <= 0) + del.u * as.numeric(x > 0)
  out <- exp(-abs(x) / coeff) / (del.l + del.u)
  return(out)
}

qasymlp <- function(p, del.l, del.u){
  cutoff <- del.l / (del.l + del.u)
  part1 <- (del.l + del.u) / (del.l * as.numeric(p <= cutoff) + del.u * as.numeric(p > cutoff))
  part2 <- ifelse(p <= cutoff, log(p), log(1 - p))
  part3 <- (del.l * as.numeric(p <= cutoff) - del.u * as.numeric(p > cutoff))
  out <- part3 * (log(part1) + part2)
  out}

rasymlp <- function(n, del.l, del.u){qasymlp(runif(n), del.l, del.u)}

#---------------
# Check the functions corresponding to AL(del.l, del.u)

del.l <- 0.3
del.u <- 0.4

R <- rasymlp(1e4, del.l, del.u)
hist(pasymlp(R, del.l, del.u), nclass = 25)

#---------------
# Fitting the full density to AL(del.l, del.u)

n <- 1e3

del.l <- 0.3
del.u <- 0.4

R <- rasymlp(n, del.l, del.u)

# maximum likelihood estimation

# true values

log(del.l) - log(1 - del.l)
log(del.u) - log(1 - del.u)

neg.log <- function(params){
  del.l <- 1 / (1 + exp(-params[1]))
  del.u <- 1 / (1 + exp(-params[2]))
  
  out <- -sum(log(dasymlp(R, del.l, del.u)))
  out}

asymlp.mles <- optim(c(-0.8472979, -0.4054651), neg.log)

del.l <- 1 / (1 + exp(-asymlp.mles$par[1]))
del.u <- 1 / (1 + exp(-asymlp.mles$par[2]))

del.l # true value = 0.3
del.u # true value = 0.4

log(0.7) - log(1 - 0.7)
log(0.9) - log(1 - 0.9)

asymlp.mles <- optim(c(0.8472979, 2.197225), neg.log)

del.l <- 1 / (1 + exp(-asymlp.mles$par[1]))
del.u <- 1 / (1 + exp(-asymlp.mles$par[2]))

del.l # true value = 0.3
del.u # true value = 0.4

# All Ok!

#---------------
# CDF and quantile function of AL(del.l, del.u) + AL(1 - del.l, 1 - del.u)

k1 <- function(del.l, del.u){
  del.l^3 / ((del.l + del.u) * (2 * del.l - 1) * (1 + del.l - del.u))}

k2 <- function(del.l, del.u){
  (del.l - 1)^3 / ((2 * del.l - 1) * (del.l - del.u - 1) * (2 - del.l - del.u))}

k3 <- function(del.l, del.u){
  del.u^3 / ((del.l + del.u) * (2 * del.u - 1) * (del.l - del.u - 1))}

k4 <- function(del.l, del.u){
  (del.u - 1)^3 / ((2 * del.u - 1) * (1 + del.l - del.u) * (2 - del.l - del.u))}

k5 <- function(del.l, del.u){2 / ((1 + 2 * del.u) * (2 * del.u - 3))}

k6 <- function(del.l, del.u){
  (-12 * del.u + 12 * del.u^2 - 5) / ((1 + 2 * del.u)^2 * (2 * del.u - 3)^2)}

k7 <- function(del.l, del.u){4 * del.u^3 / ((1 + 2 * del.u)^2 * (2 * del.u - 1))}

k8 <- function(del.l, del.u){4 * (del.u - 1)^3 / ((2 * del.u - 3)^2 * (2 * del.u - 1))}

k9 <- function(del.l, del.u){4 * del.l^3 / ((1 + 2 * del.l)^2 * (2 * del.l - 1))}

k10 <- function(del.l, del.u){4 * (del.l - 1)^3 / ((2 * del.l - 3)^2 * (2 * del.l - 1))}

k11 <- function(del.l, del.u){2 / ((1 + 2 * del.l) * (2 * del.l - 3))}

k12 <- function(del.l, del.u){
  (-12 * del.l + 12 * del.l^2 - 5) / ((1 + 2 * del.l)^2 * (2 * del.l - 3)^2)}

pasympl.sum <- function(x, del.l, del.u){
  if((del.l == 0.5) & (del.u == 0.5)){
    val <- ifelse(x <= 0, 0.5 * (1 - x) * exp(2 * x), 1 - 0.5 * (1 + x) * exp(-2 * x))
  }else if((del.l == 0.5) & (del.u != 0.5)){
    val <- ifelse(x <= 0, (k5(del.l, del.u) * x - k6(del.l, del.u)) * exp(2 * x),
                  1 - k7(del.l, del.u) * exp(-x / del.u) - k8(del.l, del.u) * exp(x / (del.u - 1)))
  }else if((del.l != 0.5) & (del.u == 0.5)){
    val <- ifelse(x <= 0, k9(del.l, del.u) * exp(x / del.l) + k10(del.l, del.u) * exp(x / (1 - del.l)),
                  1 + (k11(del.l, del.u) * x + k12(del.l, del.u)) * exp(-2 * x))
  }else{
    val <- ifelse(x <= 0, k1(del.l, del.u) * exp(x / del.l) - k2(del.l, del.u) * exp(x / (1 - del.l)),
                  1 + k3(del.l, del.u) * exp(-x / del.u) - k4(del.l, del.u) * exp(x / (del.u - 1)))
  }
  return(val)}

dasympl.sum <- function(x, del.l, del.u){
  if((del.l == 0.5) & (del.u == 0.5)){
    val <- (1 + abs(x)) * exp(-2 * abs(x)) - 0.5 * exp(-2 * abs(x))
  }else if((del.l == 0.5) & (del.u != 0.5)){
    val <- ifelse(x <= 0, exp(2 * x) * (k5(del.l, del.u) * (1 + 2 * x) - 2 * k6(del.l, del.u)),
                  k7(del.l, del.u) / del.u * exp(-x / del.u) + 
                    k8(del.l, del.u) / (1 - del.u) * exp(-x / (1 - del.u)))
  }else if((del.l != 0.5) & (del.u == 0.5)){
    val <- ifelse(x <= 0, k9(del.l, del.u) / del.l * exp(x / del.l) + 
                    k10(del.l, del.u) / (1 - del.l) * exp(x / (1 - del.l)),
                  exp(-2 * x) * (k11(del.l, del.u) * (1 - 2 * x) - 2 * k12(del.l, del.u)))
  }else{
    val <- ifelse(x <= 0, k1(del.l, del.u) / del.l * exp(x / del.l) - 
                    k2(del.l, del.u) / (1 - del.l) * exp(x / (1 - del.l)),
                  -k3(del.l, del.u) / del.u * exp(-x / del.u) + 
                    k4(del.l, del.u) / (1 - del.u) * exp(-x / (1 - del.u)))
  }
  return(val)}

qasympl.sum_0.5_0.5 <- function(p){
  solver <- function(x){abs(log(ifelse(p < 0.5, p, 1 - p)) - log(0.5) - log(1 - x) - 2 * x)}
  solver.logit <- function(xx){solver(log(xx) - log(1 - xx))}
  out <- optimize(solver.logit, c(0, 0.5), tol = 1e-10)$minimum
  out <- log(out) - log(1 - out)
  out <- ifelse(p < 0.5, out, -out)
  out}

qasympl.sum_del.l_del.u <- function(p, del.l, del.u){
  
  cdf0 <- pasympl.sum(0, del.l, del.u)
  
  if(p <= cdf0){
    solver <- function(x){
      if(del.l > 0.5){
        val <- abs(log(p) - x / del.l - 
                     log(k1(del.l, del.u) - k2(del.l, del.u) * 
                           exp(x * (2 * del.l - 1) / (del.l * (1 - del.l)))))
      }else{
        val <- abs(log(p) - x / (1 - del.l) - 
                     log(k1(del.l, del.u) * exp(x * (1 - 2 * del.l) / (del.l * (1 - del.l))) - 
                           k2(del.l, del.u)))}
      val}
    solver.logit <- function(xx){solver(log(xx) - log(1 - xx))}
    out <- optimize(solver.logit, c(0, cdf0), tol = 1e-10)$minimum
    out <- log(out) - log(1 - out)
  }else{
    solver <- function(x){
      if(del.u > 0.5){
        val <- abs(log(1 - p) + x / del.u - 
                     log(-k3(del.l, del.u) + k4(del.l, del.u) *
                           exp(x * (1 - 2 * del.u) / (del.u * (1 - del.u)))))
      }else{
        val <- abs(log(1 - p) + x / (1 - del.u) - 
                     log(-k3(del.l, del.u) * exp(x * (2 * del.u - 1) / (del.u * (1 - del.u))) + 
                           k4(del.l, del.u)))}
  val}
    solver.logit <- function(xx){solver(log(xx) - log(1 - xx))}
    out <- optimize(solver.logit, c(cdf0, 1), tol = 1e-10)$minimum
    out <- log(out) - log(1 - out)}
  out}

qasympl.sum_0.5_del.u <- function(p, del.u){
  
  del.l <- 0.5
  cdf0 <- pasympl.sum(0, del.l, del.u)
  
  if(p <= cdf0){
    solver <- function(x){
      abs(log(p) - 2 * x - log((k5(del.l, del.u) * x) - k6(del.l, del.u)))}
    solver.logit <- function(xx){solver(log(xx) - log(1 - xx))}
    out <- optimize(solver.logit, c(0, cdf0), tol = 1e-10)$minimum
    out <- log(out) - log(1 - out)
  }else{
    solver <- function(x){
      if(del.u > 0.5){
        val <- abs(log(1 - p) + x / del.u - 
                     log(k7(del.l, del.u) + k8(del.l, del.u) * 
                           exp(x * (1 - 2 * del.u) / (del.u * (1 - del.u)))))
      }else{
        val <- abs(log(1 - p) + x / (1 - del.u) - 
                     log(k7(del.l, del.u) * exp(x * (2 * del.u - 1) / (del.u * (1 - del.u))) + 
                           k8(del.l, del.u)))}
      val}
    solver.logit <- function(xx){solver(log(xx) - log(1 - xx))}
    out <- optimize(solver.logit, c(cdf0, 1), tol = 1e-10)$minimum
    out <- log(out) - log(1 - out)}
  out}

qasympl.sum_del.l_0.5 <- function(p, del.l){
  
  del.u <- 0.5
  cdf0 <- pasympl.sum(0, del.l, del.u)
  
  if(p > cdf0){
    solver <- function(x){
      abs(log(1 - p) + 2 * x - log(-k11(del.l, del.u) * x - k12(del.l, del.u)))}
    solver.logit <- function(xx){solver(log(xx) - log(1 - xx))}
    out <- optimize(solver.logit, c(cdf0, 1), tol = 1e-10)$minimum
    out <- log(out) - log(1 - out)
  }else{
    solver <- function(x){
      if(del.l > 0.5){
        val <- abs(log(p) - x / del.l - log(k9(del.l, del.u) + k10(del.l, del.u) * 
                                              exp(x * (2 * del.l - 1) / (del.l * (1 - del.l)))))
      }else{
        val <- abs(log(p) - x / (1 - del.l) - 
                     log(k9(del.l, del.u) * exp(x * (1 - 2 * del.l) / (del.l * (1 - del.l))) + 
                           k10(del.l, del.u)))}
      val}
    solver.logit <- function(xx){solver(log(xx) - log(1 - xx))}
    out <- optimize(solver.logit, c(0, cdf0), tol = 1e-10)$minimum
    out <- log(out) - log(1 - out)}
  out}

qasympl.sum <- function(p, del.l, del.u){
  if((del.l == 0.5) & (del.u == 0.5)){
    return(qasympl.sum_0.5_0.5(p))
  }else if((del.l == 0.5) & (del.u != 0.5)){
    return(qasympl.sum_0.5_del.u(p, del.u))
  }else if((del.l != 0.5) & (del.u == 0.5)){
    return(qasympl.sum_del.l_0.5(p, del.l))
  }else{
    return(qasympl.sum_del.l_del.u(p, del.l, del.u))
    }
}

#---------------
# Check the functions corresponding to AL(del.l, del.u) + AL(1 - del.l, 1 - del.u)

# Case 1: both unequal to 0.5

del.l <- 0.3
del.u <- 0.4

X <- rasymlp(1e4, del.l, del.u) + rasymlp(1e4, 1 - del.l, 1 - del.u)

hist(pasympl.sum(X, del.l, del.u), nclass = 25)

# Case 2: del.l = 0.5 only

del.l <- 0.5
del.u <- 0.3

X <- rasymlp(1e4, del.l, del.u) + rasymlp(1e4, 1 - del.l, 1 - del.u)
hist(pasympl.sum(X, del.l, del.u), nclass = 25)

# Case 3: del.u = 0.5 only

del.l <- 0.3
del.u <- 0.5

X <- rasymlp(1e4, del.l, del.u) + rasymlp(1e4, 1 - del.l, 1 - del.u)
hist(pasympl.sum(X, del.l, del.u), nclass = 25)

# Case 4: both equal to 0.5

del.l <- 0.5
del.u <- 0.5

X <- rasymlp(1e4, del.l, del.u) + rasymlp(1e4, 1 - del.l, 1 - del.u)
hist(pasympl.sum(X, del.l, del.u), nclass = 25)

#---------------
# Fitting the full density to AL(del.l, del.u) + AL(1 - del.l, 1 - del.u)

n <- 1e3

del.l <- 0.3
del.u <- 0.4

X <- rasymlp(1e4, del.l, del.u) + rasymlp(1e4, 1 - del.l, 1 - del.u)

# maximum likelihood estimation

# true values

log(del.l) - log(1 - del.l)
log(del.u) - log(1 - del.u)

neg.log <- function(params){
  del.l <- 1 / (1 + exp(-params[1]))
  del.u <- 1 / (1 + exp(-params[2]))
  
  out <- -sum(log(dasympl.sum(X, del.l, del.u)))
  out}

asymlp.sum.mles <- optim(c(-0.8472979, -0.4054651), neg.log)

del.l <- 1 / (1 + exp(-asymlp.sum.mles$par[1]))
del.u <- 1 / (1 + exp(-asymlp.sum.mles$par[2]))

del.l # true value = 0.3
del.u # true value = 0.4

log(0.7) - log(1 - 0.7)
log(0.9) - log(1 - 0.9)

neg.log <- function(params){
  del.l <- 1 / (1 + exp(-params[1]))
  del.u <- 1 / (1 + exp(-params[2]))
  
  out <- -sum(log(dasympl.sum(X, del.l, del.u)))
  out}

asymlp.sum.mles <- optim(c(0.8472979, 2.197225), neg.log)

del.l <- 1 / (1 + exp(-asymlp.sum.mles$par[1]))
del.u <- 1 / (1 + exp(-asymlp.sum.mles$par[2]))

del.l # true value = 0.3
del.u # true value = 0.4

# All OK

pasympl.sum(qasympl.sum(0.7, del.l, del.u), del.l, del.u)
pasympl.sum(qasympl.sum(0.7, del.l, del.u), 1 - del.l, 1 - del.u)
pasympl.sum(qasympl.sum(0.7, del.l, del.u), del.l, 1 - del.u)

pasympl.sum(qasympl.sum(0.7, del.l, del.u), 1 - del.l, del.u)

# Note that AL(del.l, del.u) + AL(1 - del.l, 1 - del.u) is not fully identifiable.
# AL(1 - del.l, 1 - del.u) + AL(del.l, del.u) has the same distribution.

del.l <- 0.3
del.u <- 0.4

mles <- t(sapply(1:100, function(i){
  set.seed(i)
  X <- rasymlp(1e4, del.l, del.u) + rasymlp(1e4, 1 - del.l, 1 - del.u)
  
  neg.log <- function(params){
    del.l <- 1 / (1 + exp(-params[1]))
    del.u <- 1 / (1 + exp(-params[2]))
    out <- -sum(log(dasympl.sum(X, del.l, del.u)))
    out}
  
  asymlp.sum.mles <- optim(c(1, 1), neg.log)
  del.l <- 1 / (1 + exp(-asymlp.sum.mles$par[1]))
  del.u <- 1 / (1 + exp(-asymlp.sum.mles$par[2]))
  c(del.l, del.u)}))

plot(mles)

#---------------
# bivariate PDF and CDF of X = (X1, X2)

qnorm.pasymlp <- function(x, del.l, del.u){
  coeff <- del.l * as.numeric(x <= 0) + del.u * as.numeric(x > 0)
  out <- coeff * exp(-abs(x) / coeff) / (del.l + del.u)
  out <- ifelse(x <= 0, qnorm(out), -qnorm(out))
  return(out)
}

library(mvtnorm)

pasympl.gauss <- function(x1, x2, del.l, del.u, rho){
  
  den.integrand <- function(x){
    sapply(x, function(xx){
      r <- log(xx) - log(1 - xx)
      q1 <- qnorm.pasymlp(x1 - r, 1 - del.l, 1 - del.u)
      q2 <- qnorm.pasymlp(x2 - r, 1 - del.l, 1 - del.u)
      
      q <- c(q1, q2)
      out <- pmvnorm(lower = rep(-Inf, 2), upper = q, mean = c(0, 0), 
                     corr = matrix(c(1, rho, rho, 1), nrow = 2))[1]
      out <- out * dasymlp(r, del.l, del.u)
      out <- out / (xx * (1 - xx))
      out})}
  
  out <- integrate(den.integrand, lower = 0, upper = 1)$value
  out}

dasympl.gauss <- function(x1, x2, del.l, del.u, rho){
  
  den.integrand <- function(x){
    r <- log(x) - log(1 - x)
    q1 <- qnorm.pasymlp(x1 - r, 1 - del.l, 1 - del.u)
    q2 <- qnorm.pasymlp(x2 - r, 1 - del.l, 1 - del.u)
    
    q <- cbind(q1, q2)
    out <- dmvnorm(q, mean = c(0, 0), sigma = matrix(c(1, rho, rho, 1), nrow = 2), log = T)
    out <- out + log(dasymlp(x1 - r, 1 - del.l, 1 - del.u)) +
      log(dasymlp(x2 - r, 1 - del.l, 1 - del.u)) + log(dasymlp(r, del.l, del.u))
    out <- out - dnorm(qnorm.pasymlp(x1 - r, 1 - del.l, 1 - del.u), log = T) -
                    dnorm(qnorm.pasymlp(x2 - r, 1 - del.l, 1 - del.u), log = T)
    out <- out - log(x) - log(1 - x)
    exp(out)}
  out <- integrate(den.integrand, lower = 0, upper = 1)$value
 
  out}

#---------------
# plot the density for different choices of the parameters

del.l <- 0.3
del.u <- 0.4
rho <- 0.7

x.grid <- seq(-6, 6, 0.1)
y.grid <- seq(-6, 6, 0.1)
xy.grid <- as.matrix(expand.grid(x.grid, y.grid))

dens.val <- apply(xy.grid, 1, function(x){pasympl.gauss(x[1], x[2], del.l, del.u, rho)})

library(plot3D)

scatter2D(x = xy.grid[ , 1],  y = xy.grid[ , 2], colvar = dens.val, pch = 15)

# Fitting the full density to non-unifrom-transformed data

# Simulate from a bivariate Normal copula with AL(1 - del.l, 1 - del.u) margins
# plus AL(del.l, del.u)

library(mvtnorm)

n <- 1e3

rho <- 0.2
del.l <- 0.3
del.u <- 0.4

Z <- rmvnorm(n, mean = c(0, 0), sigma = matrix(c(1, rho, rho, 1), nrow = 2))
W <- qasymlp(pnorm(Z), 1 - del.l, 1 - del.u)
R <- rasymlp(n, del.l, del.u)
X <- W + R

# check margins

hist(pasympl.sum(X[ , 1], del.l, del.u), nclass = 25)
hist(pasympl.sum(X[ , 2], del.l, del.u), nclass = 25)

# estimating del.l and del.u ignoring the joint distribution

log(del.l) - log(1 - del.l)
log(del.u) - log(1 - del.u)

neg.log <- function(params){
  del.l <- 1 / (1 + exp(-params[1]))
  del.u <- 1 / (1 + exp(-params[2]))
  
  out <- -sum(log(dasympl.sum(X[ , 1], del.l, del.u)))
  out}

asymlp.sum.mles <- optim(c(0.3, 0.2), neg.log, hessian = T)

del.l <- 1 / (1 + exp(-asymlp.sum.mles$par[1]))
del.u <- 1 / (1 + exp(-asymlp.sum.mles$par[2]))

del.l # true value = 0.3/0.7
del.u # true value = 0.4/0.6

# maximum likelihood estimation for joint

transform <- list(
  logit = function(x, lower = 0, upper = 1){
    x <- (x - lower) / (upper - lower)
    return(log(x / (1 - x)))
  },
  inv.logit = function(x, lower = 0, upper = 1){
    p <- exp(x) / (1 + exp(x))
    p <- p * (upper - lower) + lower
    return(p)
  }
)

# true values

n <- 5e2

rho <- 0.2
del.l <- 0.3
del.u <- 0.4

Z <- rmvnorm(n, mean = c(0, 0), sigma = matrix(c(1, rho, rho, 1), nrow = 2))
W <- qasymlp(pnorm(Z), 1 - del.l, 1 - del.u)
R <- rasymlp(n, del.l, del.u)
X <- W + R

init <- c(transform$logit(del.l, 0.01, 0.99), 
          transform$logit(del.u, 0.01, 0.99), 
          transform$logit(rho, -0.99, 0.99))

neg.log <- function(params){
  del.l <- transform$inv.logit(params[1], 0.01, 0.99)
  del.u <- transform$inv.logit(params[2], 0.01, 0.99)
  rho <- transform$inv.logit(params[3], -0.99, 0.99)
  
  out <- sum(sapply(1:nrow(X), function(i){
    -log(dasympl.gauss(X[i, 1], X[i, 2], del.l, del.u, rho))}))
  out}

mles <- optim(init, neg.log, hessian = T)

transform$inv.logit(mles$par[1], 0.01, 0.99)
transform$inv.logit(mles$par[2], 0.01, 0.99)
transform$inv.logit(mles$par[3], -0.99, 0.99)

#--------------------
# copula part

library(mvtnorm)

n <- 1e2

rho <- 0.2
del.l <- 0.3
del.u <- 0.4

Z <- rmvnorm(n, mean = c(0, 0), sigma = matrix(c(1, rho, rho, 1), nrow = 2))
W <- qasymlp(pnorm(Z), 1 - del.l, 1 - del.u)
R <- rasymlp(n, del.l, del.u)
X <- W + R
X <- apply(X, 2, rank) / (nrow(X) + 1)

c.u1.u2 <- function(u1, u2, del.l, del.u, rho){
  x.1 <- qasympl.sum(u1, del.l, del.u)
  x.2 <- qasympl.sum(u2, del.l, del.u)
  out <- dasympl.gauss(x.1, x.2, del.l, del.u, rho) / 
    (dasympl.sum(x.1, del.l, del.u) * dasympl.sum(x.2, del.l, del.u))
  out}

# plot the copula values across the unit square

del.l <- 0.3
del.u <- 0.4
rho <- 0.7

u1.grid <- seq(0.01, 0.99, 0.01)
u2.grid <- seq(0.01, 0.99, 0.01)
u.grid <- as.matrix(expand.grid(u1.grid, u2.grid))

copula.dens.val <- apply(u.grid, 1, function(u){c.u1.u2(u[1], u[2], del.l, del.u, rho)})

library(plot3D)

scatter2D(x = u.grid[ , 1],  y = u.grid[ , 2], colvar = log(copula.dens.val), pch = 15)

# log-likelihood as a function of each parameter

n <- 1e2

rho <- 0.2
del.l <- 0.3
del.u <- 0.4

Z <- rmvnorm(n, mean = c(0, 0), sigma = matrix(c(1, rho, rho, 1), nrow = 2))
W <- qasymlp(pnorm(Z), 1 - del.l, 1 - del.u)
R <- rasymlp(n, del.l, del.u)
X <- W + R
X <- apply(X, 2, rank) / (nrow(X) + 1)

ll.del.l <- function(del.l){
  del.u <- 0.4
  rho <- 0.2
  
  out <- sum(sapply(1:nrow(X), function(i){
    log(c.u1.u2(X[i, 1], X[i, 2], del.l, del.u, rho))}))
  out}

ll.vals <- sapply(seq(0.2, 0.6, 0.01), ll.del.l)

plot(seq(0.2, 0.6, 0.01), exp(ll.vals), type = "l")

ll.del.u <- function(del.u){
  del.l <- 0.3
  rho <- 0.2
  
  out <- sum(sapply(1:nrow(X), function(i){
    log(c.u1.u2(X[i, 1], X[i, 2], del.l, del.u, rho))}))
  out}

ll.vals <- sapply(seq(0.2, 0.6, 0.01), ll.del.u)

plot(seq(0.2, 0.6, 0.01), ll.vals, type = "l")

#3. negative likelihood

n <- 1e2

rho <- 0.2
del.l <- 0.7
del.u <- 0.4

Z <- rmvnorm(n, mean = c(0, 0), sigma = matrix(c(1, rho, rho, 1), nrow = 2))
W <- qasymlp(pnorm(Z), 1 - del.l, 1 - del.u)
R <- rasymlp(n, del.l, del.u)
X <- W + R
X <- apply(X, 2, rank) / (nrow(X) + 1)
#X <- cbind(pasympl.sum(X[ , 1], del.l, del.u), pasympl.sum(X[ , 2], del.l, del.u))


init <- c(transform$logit(del.l, 0.01, 0.99), 
          transform$logit(del.u, 0.01, 0.99), 
          transform$logit(rho, -0.99, 0.99))

neg.log <- function(params){
 del.l <- transform$inv.logit(params[1], 0.01, 0.99)
  del.u <- transform$inv.logit(params[2], 0.01, 0.99)
  rho <- transform$inv.logit(params[3], -0.99, 0.99)
  
  out <- sum(sapply(1:nrow(X), function(i){
    -log(c.u1.u2(X[i, 1], X[i, 2], del.l, del.u, rho))}))
  out}

mles <- optim(init, neg.log, hessian = T)
#mles <- optim(c(0.1,0.2,0.3), neg.log, hessian = T)

transform$inv.logit(mles$par[1], 0.01, 0.99)
transform$inv.logit(mles$par[2], 0.01, 0.99)
transform$inv.logit(mles$par[3], -0.99, 0.99)

# sqrt(diag(solve(mles$hessian)))

sqrt(diag(solve(mles$hessian))[1] / transform$inv.logit(mles$par[1], 0.01, 0.99)^2)
sqrt(diag(solve(mles$hessian))[2] / transform$inv.logit(mles$par[2], 0.01, 0.99)^2)
sqrt(diag(solve(mles$hessian))[3] / transform$inv.logit(mles$par[3], 0.01, 0.99)^2) 
#############################################################################################
#install.packages("ismev")
#install.packages("evd")
#install.packages("readr")


library(ismev)
library(readr)
library(evd)
library(mvtnorm)

Indigo <- read_csv("INGL Historical Data.csv")
SpiceJet <- read_csv("SPJT Historical Data.csv")

Indigo <- Indigo[-c(497:507), ]
SpiceJet <- SpiceJet[-c(498:508), ]

years.Indigo <- as.vector(sapply(Indigo$Date, function(x){
  as.numeric(strsplit(x, split = " ")[[1]][3])}))

years.SpiceJet <- as.vector(sapply(SpiceJet$Date, function(x){
  as.numeric(strsplit(x, split = " ")[[1]][3])}))

days.Indigo <- as.vector(sapply(Indigo$Date, function(x){
  as.numeric(strsplit(strsplit(x, split = ", ")[[1]][1], split = " ")[[1]][2])}))

days.SpiceJet <- as.vector(sapply(SpiceJet$Date, function(x){
  as.numeric(strsplit(strsplit(x, split = ", ")[[1]][1], split = " ")[[1]][2])}))

months.Indigo <- as.vector(sapply(Indigo$Date, function(x){strsplit(x, split = " ")[[1]][1]}))

months.SpiceJet <- as.vector(sapply(SpiceJet$Date, function(x){strsplit(x, split = " ")[[1]][1]}))

months.Indigo <- as.vector(sapply(months.Indigo, function(x){which(month.abb == x)}))

months.SpiceJet <- as.vector(sapply(months.SpiceJet, function(x){which(month.abb == x)}))

identify.issue <- which(sapply(1:497, function(i){
  sum((months.Indigo - months.SpiceJet[-i])^2) +
    sum((years.Indigo - years.SpiceJet[-i])^2) +
    sum((days.Indigo - days.SpiceJet[-i])^2)}) == 0)

SpiceJet <- SpiceJet[-identify.issue, ]






# data preprocessing 

n <- dim(Indigo)[1]
likelihood.before <- numeric(length = n-20 +1)
likelihood.after <-  numeric(length = n-20 +1)
mles.before.store <-  matrix(0,n-19,3)
mles.after.store <-  matrix(0,n-19,3)
n<- dim(Indigo)[1] 
i =10
tick <- proc.time()[3]
for(i in 10:(n-10)){
  ingl.split = split(Indigo, f = (1:n) >= i)
  spcjt.split = split(SpiceJet, f = (1:n) >= i)
  ingl.before = ingl.split$'FALSE'$Price
  ingl.after = ingl.split$'TRUE'$Price
  spcjt.before = spcjt.split$'FALSE'$Price
  spcjt.after = spcjt.split$'TRUE'$Price
  X <- cbind(ingl.before, spcjt.before)
  #X <- apply(X, 2, rank) / (nrow(X) + 1) 
  emp.para <- gev.fit(X[,1])$mle 
  X[,1] <- pgev(X[,1], loc =  emp.para[1], scale = emp.para[2] , shape = emp.para[3] )
  emp.para <- gev.fit(X[,2])$mle 
  X[,2] <- pgev(X[,2], loc =  emp.para[1], scale = emp.para[2] , shape = emp.para[3] )
  #############remove values very near to 1 and 0 ##############
  ########### chosen 1e-8 since 1-1e-8 is nearly considered to be in r
  identify.issue <- which(1-X[,2] < 1e-8)
  if(length(identify.issue) > 0){
    X <- X[-identify.issue ,]}
  identify.issue <- which(1-X[,1] < 1e-8)
  if(length(identify.issue) > 0){
    X <- X[-identify.issue ,]}
  
  ############################################################
  init <- c(transform$logit(0.2, 0.01, 0.99), 
            transform$logit(0.2, 0.01, 0.99), 
            transform$logit(0.2, -0.99, 0.99))
  mles.before <- optim(init, neg.log, hessian = T)$par
  mles.before <- c(transform$inv.logit(mles.before[1], 0.01, 0.99),
                   transform$inv.logit(mles.before[2], 0.01, 0.99),
                   transform$inv.logit(mles.before[3], -0.99, 0.99))
  mles.before.store[i,] <- mles.before
  likelihood.before[i-10 +1] <- sum(sapply(1:nrow(X), function(i){
    log(c.u1.u2(X[i, 1], X[i, 2], mles.before[1],mles.before[2] , mles.before[3]))}))
  
  
  
  init <- c(transform$logit(0.3, 0.01, 0.99), 
            transform$logit(0.4, 0.01, 0.99), 
            transform$logit(0.2, -0.99, 0.99))
  X <- cbind(ingl.after, spcjt.after)

  X <- apply(X, 2, rank) / (nrow(X) + 1) 
#  emp.para <- gev.fit(X[,1])$mle 
#  X[,1] <- pgev(X[,1], loc =  emp.para[1], scale = emp.para[2] , shape = emp.para[3] )
#  emp.para <- gev.fit(X[,2])$mle 
#  X[,2] <- pgev(X[,2], loc =  emp.para[1], scale = emp.para[2] , shape = emp.para[3] )
  #############remove values very near to 1 and 0 ##############
  ########### chosen 1e-8 since 1-1e-8 is nearly considered to be in r
  identify.issue <- which(1-X[,2] < 1e-8)
  if(length(identify.issue) > 0){
    X <- X[-identify.issue ,]}
  identify.issue <- which(1-X[,1] < 1e-8)
  if(length(identify.issue) > 0){
    X <- X[-identify.issue ,]}
  
  ############################################################
  mles.after <- optim(init, neg.log, hessian = T)$par 
  mles.after <- c(transform$inv.logit(mles.after[1], 0.01, 0.99),
                   transform$inv.logit(mles.after[2], 0.01, 0.99),
                   transform$inv.logit(mles.after[3], -0.99, 0.99))
  mles.after.store[i,] <- mles.after
  likelihood.after[i-10 +1] <- sum(sapply(1:nrow(X), function(i){
    log(c.u1.u2(X[i, 1], X[i, 2], mles.after[1],mles.after[2] , mles.after[3]))}))
} 

tock <- proc.time()[3]

minutes = (tock - tick) / 60



#############################################33#whole data
init <- c(transform$logit(0.2, 0.01, 0.99), 
          transform$logit(0.2, 0.01, 0.99), 
          transform$logit(0.2, -0.99, 0.99))
X <- cbind(Indigo$Price, SpiceJet$Price)
# X <- apply(X, 2, rank) / (nrow(X) + 1) 
emp.para <- gev.fit(X[,1])$mle 
X[,1] <- pgev(X[,1], loc =  emp.para[1], scale = emp.para[2] , shape = emp.para[3] )
emp.para <- gev.fit(X[,2])$mle 
X[,2] <- pgev(X[,2], loc =  emp.para[1], scale = emp.para[2] , shape = emp.para[3] )
 

init <- c(transform$logit(0.2, 0.01, 0.99), 
          transform$logit(0.2, 0.01, 0.99), 
          transform$logit(0.2, -0.99, 0.99))
mles <- optim(init, neg.log, hessian = T)$par 
mles <- c(transform$inv.logit(mles$par[1], 0.01, 0.99),
                transform$inv.logit(mles$par[2], 0.01, 0.99),
                transform$inv.logit(mles$par[3], -0.99, 0.99))

likelihood.whole<- sum(sapply(1:nrow(X), function(i){
  log(c.u1.u2(X[i, 1], X[i, 2], mles[1],mles[2] , mles[3]))}))

change.ratio <- (likelihood.after*likelihood.before)/likelihood.whole
changepoint <- Indigo[which.max(change.ratio) + 9]$ï..Date
plot(change.ratio)



###########################3simulation study

data.lengths = c(1000) 
reps = 10

para.store1<-matrix(rep(list(), 3),nrow = reps, ncol = length(data.lengths)) 

minutes1 <- numeric(length = reps)

for(i in 1:reps){
  tick <- proc.time()[3]
  for(j in 1:length(data.lengths)){ 
    set.seed(i)
n <- data.lengths[j]

rho <- 0.2
del.l <- 0.7
del.u <- 0.8

Z <- rmvnorm(n, mean = c(0, 0), sigma = matrix(c(1, rho, rho, 1), nrow = 2))
W <- qasymlp(pnorm(Z), 1 - del.l, 1 - del.u)
R <- rasymlp(n, del.l, del.u)
X <- W + R
X <- apply(X, 2, rank) / (nrow(X) + 1)
#X <- cbind(pasympl.sum(X[ , 1], del.l, del.u), pasympl.sum(X[ , 2], del.l, del.u))


init <- c(transform$logit(0.5, 0.01, 0.99), 
          transform$logit(0.5, 0.01, 0.99), 
          transform$logit(0.5, -0.99, 0.99))

mles <- optim(init, neg.log, hessian = T)

para.store1[[reps*(j-1) + i]] <- c(transform$inv.logit(mles$par[1], 0.01, 0.99),
                    transform$inv.logit(mles$par[2], 0.01, 0.99),
                   transform$inv.logit(mles$par[3], -0.99, 0.99))
  }
  tock <- proc.time()[3]
  minutes1[i] = (tock - tick) / 60
}

save.image(file = "my_work_space020708.RData")
save.image(file = "my1000_work_space020708.RData")

