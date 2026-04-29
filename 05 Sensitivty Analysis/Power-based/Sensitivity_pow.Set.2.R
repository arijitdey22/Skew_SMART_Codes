rm(list = ls())

library(sn)
library(fGarch)

###------------------------------
###calculating sample sizes

samp.size.pow.s2 <- function(gam.11 = 0.35, gam.12 = 0.35, nu.d1 = 1, eff.size = 0.2,
                              pow = 0.8, alpha = 0.05, r1 = 1, r2 = 1)
{
  xi.d1 <- nu.d1 / sqrt(1+nu.d1^2)
  kap.d1 <- (1 - 2*xi.d1^2 / pi) 
  ret <- 2 * (qnorm(pow) + qnorm(1- alpha/2))^2 * ((4 - 2*gam.11) + r1 * r2 * (4 - 2*gam.12)) * kap.d1 / (eff.size^2 * (1+r1))
  return(ceiling(ret))
}

#===============================================================================
#drawing random samples from truncated skew-normal distribution

r.trun.sn <- function(n, mu = 0, sigma = 1, lambda = 0, lower = -Inf, upper = Inf) {
  
  out <- numeric(n)
  oversample_factor <- 1 / (psn(upper, mu, sigma, lambda) - psn(lower, mu, sigma, lambda))
  
  i <- 1
  while (i <= n) {
    m <- ceiling((n - i + 1) * oversample_factor * 1.2)
    x <- rsn(m, mu, sigma, lambda)
    x <- x[x > lower & x < upper]
    n_add <- min(length(x), n - i + 1)
    if (n_add > 0) {
      out[i:(i + n_add - 1)] <- x[1:n_add]
      i <- i + n_add
    }
  }
  return(out)
}

###------------------------------
###calculating nu-relationship from kappa-relationship
#input: nu_1; output: nu_2

nu.2 <- function(nu.1, r){
  
  sign.nu <- ifelse(nu.1 >= 0, 1, -1)
  A = pi * (1 - r) / 2 + nu.1^2 / (1+nu.1^2)
  nu.2 = sign.nu * sqrt(abs(A / (1 - A)))
  
  return(nu.2)
  
}

#===============================================================================
#the main function

sim.fun.pow.set.2 <- function(gam.11, gam.12, nu.d1, nu.samp, eff.size, pow = 0.8, alpha = 0.05,
                         r1 = 1, r2 = 1, rep = 1e4, N = NA, pi.11 = 0.5, pi.22 = 0.5)
{
  
  count <- 0                              
  pow.cap.del.a <- numeric()                    #.a denotes 'all'
  eff.size.cap.a <- numeric()
  
  #1------
  
  if (is.na(N)) {N = samp.size.pow.s2(gam.11, gam.12, nu.samp, eff.size, pow)}
  
  #-------
  
  for (i in 1:rep)
  {
    set.seed(i)
    
    #2------
    mu.l.1 <- runif(1, 0.1, 5)
    sigma.l.1 <- runif(1, 0.1,5)
    nu.l.1 <- runif(1, 0.1,5)
    
    ita <- mu.l.1 + sigma.l.1 * qsn((1 - gam.11), 0, 1, nu.l.1)
    
    sigma.l.2 <- sqrt(r1) * sigma.l.1
    nu.l.2 <- nu.2(nu.l.1, r2)
    mu.l.2 <- ita - sigma.l.2 * qsn((1 - gam.12), 0, 1, nu.l.2)
    
    #3------
    N.1 <- ceiling(N * pi.11)
    if (N.1 == 0){next}
    
    L.1 <- as.numeric(rsn(N.1, mu.l.1, sigma.l.1, nu.l.1))
    
    # ~ ~ ~ ~ ~
    
    N.2 <- ceiling(N * pi.11)
    if (N.2 == 0) {next}
    
    L.2 <- as.numeric(rsn(N.2, mu.l.2, sigma.l.2, nu.l.2))
    
    #4------
    N.1.R <- sum(L.1 > ita)
    N.1.NR <- N.1 - N.1.R
    if (N.1.R == 0 | N.1.NR == 0) {next}
    
    # ~ ~ ~ ~ ~
    
    N.2.R <- sum(L.2 > ita)
    N.2.NR <- N.2 - N.2.R
    if (N.2.R == 0 | N.2.NR == 0) {next}
    
    #5------
    m.1.R <- r.trun.sn(N.1.R, mu.l.1, sigma.l.1, nu.l.1, ita, Inf)
    m.1.NR <- r.trun.sn(N.1.NR, mu.l.1, sigma.l.1, nu.l.1, -Inf, ita)
    
    if (length(m.1.R)  != N.1.R | length(m.1.NR)  != N.1.NR) {next}
    foo <- snormFit(m.1.R)$par
    mu.m.1.R <- foo[[1]]
    
    foo <- snormFit(m.1.NR)$par
    mu.m.1.NR <- foo[[1]]
    
    #6------
    N.1.NR.2 <- ceiling(pi.22 * N.1.NR)
    if (N.1.NR.2 == 0) {next}
    
    # ~ ~ ~ ~ ~
    
    N.2.NR.2 <- ceiling(pi.22 * N.2.NR)
    if (N.2.NR.2 == 0) {next}
    
    #7------
    mu.1.R.1 <- 0.1 + log(2) * mu.l.1 + 0.1 * mu.m.1.R
    mu.1.NR.2 <- 0.1 + log(2) * mu.l.1 + 0.5 * mu.m.1.NR
    
    # ~ ~ ~ ~ ~
    # later
    
    #8------
    mu.d1 <- gam.11 * mu.1.R.1 + (1-gam.11) * mu.1.NR.2
    
    # ~ ~ ~ ~ ~
    sigma.d1 <- runif(1, 3, 6)
    sigma.d3 <- sqrt(r1) * sigma.d1
    
    mu.d3 <- mu.d1 - eff.size * sqrt((sigma.d1^2 + sigma.d3^2) / 2)
    #mu.d3 <- mu.d1 - eff.size
    
    mu.2.R.1 <- mu.d3 + (mu.1.R.1 - mu.d1)
    mu.2.NR.2 <- (mu.d3 - gam.12 * mu.2.R.1) / (1-gam.12)
    
    #9------
    
    Y.1.R.1 <- as.numeric(rsn(N.1.R, mu.1.R.1, sigma.d1, nu.d1))
    Y.1.NR.2 <- as.numeric(rsn(N.1.NR.2, mu.1.NR.2, sigma.d1, nu.d1))
    
    # ~ ~ ~ ~ ~
    
    nu.d3 <- nu.2(nu.d1, r2)
    
    Y.2.R.1 <- as.numeric(rsn(N.2.R, mu.2.R.1, sigma.d3, nu.d3))
    Y.2.NR.2 <- as.numeric(rsn(N.2.NR.2, mu.2.NR.2, sigma.d3, nu.d3))
    
    #10------
    Y.bar.d1 <- (sum(2* Y.1.R.1) + sum(4 * Y.1.NR.2)) / (2*length(Y.1.R.1) + 4*length(Y.1.NR.2))
    var.cap.Y.bar.d1 <- (   sum(4 * (Y.1.R.1 - Y.bar.d1)^2) + sum(16 * (Y.1.NR.2 - Y.bar.d1)^2)   ) / N^2
    
    # ~ ~ ~ ~ ~
    
    Y.bar.d3 <- (sum(2* Y.2.R.1) + sum(4 * Y.2.NR.2)) / (2*length(Y.2.R.1) + 4*length(Y.2.NR.2))
    var.cap.Y.bar.d3 <- (   sum(4 * (Y.2.R.1 - Y.bar.d3)^2) + sum(16 * (Y.2.NR.2 - Y.bar.d3)^2)   ) / N^2
    
    del.cap <- (Y.bar.d1 - Y.bar.d3) - ( sigma.d1 * ( nu.d1 / sqrt(1 + nu.d1^2) ) - sigma.d3 * ( nu.d3 / sqrt(1 + nu.d3^2) ) ) * sqrt(2/pi)
    
    #11------
    Z.2.cap <- del.cap / sqrt(var.cap.Y.bar.d1 + var.cap.Y.bar.d3)
    
    #12------
    pow.cap.del <- (Z.2.cap < - qnorm(1 -alpha / 2)) | ((Z.2.cap > qnorm(1- alpha / 2)))
    
    #---------------------
    #store
    pow.cap.del.a <- c(pow.cap.del.a, pow.cap.del)
    eff.size.cap.a <- c(eff.size.cap.a, del.cap / (sqrt(sigma.d1^2 + sigma.d3^2/2)) )
    
    #---------------------
    #aesthetics
    if ((i %% 100 == 0) & (i < rep))
    {
      print(paste0("We just finished ", i, " steps."))
    }
    
    count <- count + 1
    
  }
  
  #------------------------------------
  
  print(paste0("We just finished ", rep, " steps."))
  
  ret <- list(Power = mean(na.omit(pow.cap.del.a)), Eff.size_hat = mean(na.omit(eff.size.cap.a)),
              Sample_Size = N, Effective_loops = count)
  
  return(ret)
  
}

#===============================================================================
#tabulating values

eff.size <- 0.15
gam.11 <- 0.5
gam.12 <- 0.55

nu.vec <- c(1, 5)
c.nu.vec <- c(0.75, 0.9, 1.1, 1.25)

store.pow <- matrix(0, ncol = length(c.nu.vec), nrow = length(nu.vec))

for (i.nu in 1:length(nu.vec)){
  for (i.c.nu in 1:length(c.nu.vec)){
    nu.true <- nu.vec[i.nu]
    nu.samp <- nu.true * c.nu.vec[i.c.nu]
    foo <- sim.fun.pow.set.2(gam.11, gam.12, nu.true, nu.samp, eff.size)
    
    store.pow[i.nu, i.c.nu] <- foo$Power
  }
}

store.pow <- as.table(store.pow)
colnames(store.pow) <- c("c.nu = 0.75", "c.nu = 0.9", "c.nu = 1.1", "c.nu = 1.25")
rownames(store.pow) <- c("nu = 1", "nu = 5")
(store.pow <- round(store.pow, 3))

save(store.pow, file = "Sensitivity_pow.Set.2.Rdata")
