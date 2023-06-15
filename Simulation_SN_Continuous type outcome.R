library(stats4)
library(sn)                   #for skew normal
library(moments)              #for skewness of a sample

######=========================================================================================================
######=========================================================================================================

#rm(list = ls())

samp.size.set1.sn <- function(gam.a, psi, nu.a, alpha = 0.05)
{
  xi.a <- nu.a / sqrt(1+nu.a^2)
  ret <- (4 - 2*gam.a) * qnorm(alpha/2)^2 * (1 - 2*xi.a^2 / pi) / psi^2
  return(ceiling(ret))
}

r.trun.sn <- function(n, mu, sigma, lambda, lower, upper)   #requires sn package
{
  samp <- numeric()
  for (i in 1:n)
  {
    accept <- 0
    try <- 0
    while ((accept == 0) & (try < 1e3))
    {
      y <- rsn(1, mu, sigma, lambda)
      
      if ((y < upper) & (y > lower))
      {
        accept <- 1
        samp <- c(samp, y)
      }
      try <- try + 1
    }
  }
  
  return(samp)
}

par.estimate.sn <- function(vec)
{
  skw <- min(0.99, abs(skewness(vec)))
  del <- sqrt((pi / 2) * skw^(2/3) / (  skw^(2/3) + (2 - pi/2)^(2/3))  )
  sign <- ifelse(skw > 0, 1, -1)
  nu <- sign * del / sqrt(1- del^2)            #sometimes del.1 is coming to be greater than 1
  
  #----
  
  del <- nu / sqrt(1 + nu^2)
  sigma.sq <- var(vec) / (1 - 2*del^2/pi)
  sigma <- sqrt(sigma.sq)
  
  #----
  
  mu <- mean(vec) - sigma * del * sqrt(2/pi)
  
  #----

  ret <- list(mu, sigma, nu)
  return(ret)
}


cont.set1.sn <- function(psi, gam.a, nu.a, rep = 3e3, alpha = 0.05, pi.a = 0.5, pi.f = 0.5)
{
  
  count <- 0                                #count for effective loops (that are not terminated by 'next')
  CP.del1.a <- numeric()                    #.a represents all i.e. all the reps
  psi.cap.del1.a <- numeric()
  
  #1
  #N = samp.size.set1.sn(gam.a, psi, nu.a)
  N = 1e3
  
  for (i in 1:rep)
  {
    
    #2
    mu.l.a <- runif(1, 0.1, 5)
    sigma.l.a <- runif(1, 0.1,5)
    #~~~~~~~~~~~~~~~~~~
    nu.l.a <- runif(1, 0.1,5)
    #~~~~~~~~~~~~~~~~~~
    
    ita <- mu.l.a + sigma.l.a * qsn((1 - gam.a), 0, 1, nu.l.a)
    
    #3
    N.a <- ceiling(N*pi.a)
    if (N.a == 0)
    {
      next
    }
    
    samp1 <- as.numeric(rsn(N.a, mu.l.a, sigma.l.a, nu.l.a))
    
    foo <- par.estimate.sn(samp1)
    
    mu.l.a <- foo[[1]]
    sigma.l.a <- foo[[2]]
    nu.l.a <- foo[[3]]
   
    #4
    N.a.R <- sum(samp1 > ita)
    if (N.a.R == 0)
    {
      next
    }
    N.a.NR <- N.a - N.a.R
    if (N.a.NR == 0)
    {
      next
    }

    #5
    L2.a.R <- r.trun.sn(N.a.R, mu.l.a, sigma.l.a, nu.l.a, ita, Inf)
    L2.a.NR <- r.trun.sn(N.a.NR, mu.l.a, sigma.l.a, nu.l.a, -Inf, ita)
    
    if(length(L2.a.R)  != N.a.R )
    {
      next
    }
    if (length(L2.a.NR)  != N.a.NR)
    {
      next
    }
    
    foo <- par.estimate.sn(L2.a.R)
    mu.l.a.R <- foo[[1]]
    sigma.l.a.R <- foo[[2]]
    nu.l.a.R <- foo[[3]]
    
    foo <- par.estimate.sn(L2.a.NR)
    mu.l.a.NR <- foo[[1]]
    sigma.l.a.NR <- foo[[2]]
    nu.l.a.NR <- foo[[3]]
    
    #6
    N.a.NR.f <- ceiling(pi.f * N.a.NR)
    if (N.a.NR.f == 0)
    {
      next
    }
  
    #7
    mu.a.R <- 0.1 + log(2) * mu.l.a + 0.1 * mu.l.a.R
    mu.a.NR.f <- 0.1 + log(2) * mu.l.a + 0.5 * mu.l.a.NR
    
    #8
    mu.d1 <- gam.a * mu.a.R + (1-gam.a)*mu.a.NR.f
    
    del1 <- mu.d1
    
    #9
    sigma.d1 <- runif(1, 1, 10)
    #~~~~~~~~~~~~~~~~~~
    nu.d1 <- runif(1, 1, 10)
    #~~~~~~~~~~~~~~~~~~
    
    Y.a.R <- as.numeric(rsn(N.a.R, mu.a.R, sigma.d1, nu.d1))
    Y.a.NR.f <- as.numeric(rsn(N.a.NR.f, mu.a.NR.f, sigma.d1, nu.d1))
    
    #10
    Y.bar.d1 <- (sum(2* Y.a.R) + sum(4 * Y.a.NR.f)) / (2*length(Y.a.R) + 4*length(Y.a.NR.f))
    #~~~~~~~~~~~~~~~~~~
    var.cap.Y.bar.d1 <- (sum(4 * (Y.a.R - Y.bar.d1)^2) + sum(16 * (Y.a.NR.f - Y.bar.d1)^2))  / N^2
    #xi.d1 <- nu.d1 / sqrt(1+nu.d1^2)
    #var.cap.Y.bar.d1 <- ((sum(4 * (Y.a.R - Y.bar.d1)^2) + sum(16 * (Y.a.NR.f - Y.bar.d1)^2))  / N^2) * (1 - 2*xi.d1^2 / pi)
    #~~~~~~~~~~~~~~~~~~
    
    #11
    del1.cap <- Y.bar.d1
    w1.cap <- qnorm(1 - alpha / 2) * sqrt(var.cap.Y.bar.d1)
    
    #---------------------
    
    CP.del1.a <- c(CP.del1.a, ((del1 <= del1.cap + w1.cap) & (del1 >= del1.cap - w1.cap)))
    psi.cap.del1.a <- c(psi.cap.del1.a, (w1.cap / sigma.d1))
    
    if ((i %% 1000 == 0) & (i < rep))
    {
      print(paste0("We just finished ", i/1000, "k steps."))
    }
    
    count <- count + 1
  }
  
  print(paste0("We just finished ", rep/1000, "k steps."))
  ret <- list(Coverage_Prob = mean(CP.del1.a), pis_hat = mean(psi.cap.del1.a), Effective_loops = count)
  
  return(ret)
}


#cont.set1.sn(psi, gam.a, nu.a)

cont.set1.sn(0.35, 0.35, 0.35)
cont.set1.sn(0.5, 0.35, 0.35)
cont.set1.sn(0.8, 0.35, 0.35)

cont.set1.sn(0.35, 0.35, 0.65)
cont.set1.sn(0.5, 0.35, 0.65)
cont.set1.sn(0.8, 0.35, 0.65)


a


























###=========================================================================================================
###=========================================================================================================
###=========================================================================================================

######=========================================================================================================
######=========================================================================================================

#rm(list = ls())

samp.size.set2.sn <- function(gam.a, psi, nu.a, alpha = 0.05)
{
  #need to calculate
}

r.trun.sn <- function(n, mu, sigma, lambda, lower, upper)   #requires sn package
{
  samp <- numeric()
  for (i in 1:n)
  {
    accept <- 0
    try <- 0
    while ((accept == 0) & (try < 1e3))
    {
      y <- rsn(1, mu, sigma, lambda)
      
      if ((y < upper) & (y > lower))
      {
        accept <- 1
        samp <- c(samp, y)
      }
      try <- try + 1
    }
  }
  
  return(samp)
}

nu.estimate <- function(vec)
{
  skw <- min(0.99, abs(skewness(vec)))
  del <- sqrt((pi / 2) * skw^(2/3) / (skw^(2/3) + (2 - pi/2)^(2/3)))
  sign <- ifelse(skw > 0, 1, -1)
  nu <- sign * del / sqrt(1- del^2)            #sometimes del.1 is coming to be greater than 1
  return(nu)
}

cont.set1.sn












for (i in 1:rep)
{
  
  #1
  r.1 = 1
  r.2 = 1
  alpha = 0.05
  gam.a = gam.b = 0.35
  pi.a = 0.5
  pi.f = 0.25
  #N = samp.size.set2(0.35, 0.5, 1)
  N = sample(100:150, 1)
  
  #========================================================================
  
  #2
  mu.l.a <- runif(1, 0.1, 5)
  sigma.l.a <- runif(1, 0.1,5)
  nu.l.a <- runif(1, 0.1, 5)
  
  ita <- mu.l.a + sigma.l.a * qsn(1 - gam.a, 0, 1, nu.l.a)
  
  sigma.l.b <- sqrt(r.1 * sigma.l.a^2)
  nu.l.b <- r.2 * nu.l.a
  
  mu.l.b <- ita - sigma.l.b * qsn(1 - gam.b, 0, 1, nu.l.b)
  
  #========================================================================
  
  #3
  N.a <- ceiling(N*pi.a)
  N.b <- N - N.a
  samp1 <- as.numeric(rsn(N.a, mu.l.a, sigma.l.a, nu.l.a))
  samp2 <- as.numeric(rsn(N.b, mu.l.b, sigma.l.b, nu.l.b))
  
  #new mu's
  mu.l.a <- mean(samp1)
  mu.l.b <- mean(samp2)
  
  #new sigma's
  sigma.l.a <- sqrt(var(samp1))
  sigma.l.b <- sqrt(var(samp2))
  
  #new nu's
  nu.l.a <- nu.estimate(samp1)
  nu.l.b <- nu.estimate(samp2)
  
  #========================================================================
  
  #4
  N.a.R <- sum(samp1 > ita)
  N.a.NR <- N.a - N.a.R
  
  N.b.R <- sum(samp2 > ita)
  N.b.NR <- N.b - N.b.R
  
  #========================================================================
  
  #5
  L2.a.R <- r.trun.sn(N.a.R, mu.l.a, sigma.l.a, nu.l.a, ita, Inf)
  L2.a.NR <- r.trun.sn(N.a.NR, mu.l.a, sigma.l.a, nu.l.a, -Inf, ita)
  L2.b.R <- r.trun.sn(N.b.R, mu.l.b, sigma.l.b, nu.l.b, ita, Inf)
  L2.b.NR <- r.trun.sn(N.b.NR, mu.l.b, sigma.l.b, nu.l.b, -Inf, ita)

  mu.l.a.R <- mean(L2.a.R)
  sigma.l.a.R <- sqrt(var(L2.a.R))
  nu.l.a.R <- nu.estimate(L2.a.R)
  
  mu.l.a.NR <- mean(L2.a.NR)
  sigma.l.a.NR <- sqrt(var(L2.a.NR))
  nu.l.a.NR <- nu.estimate(L2.a.NR)
  
  mu.l.b.R <- mean(L2.b.R)
  sigma.l.b.R <- sqrt(var(L2.b.R))
  nu.l.b.R <- nu.estimate(L2.b.R)
  
  mu.l.b.NR <- mean(L2.b.NR)
  sigma.l.b.NR <- sqrt(var(L2.b.NR))
  nu.l.b.NR <- nu.estimate(L2.b.NR)
  
  #========================================================================
  
  #6
  N.a.NR.f <- ceiling(pi.f * N.a.NR)
  N.b.NR.f <- ceiling(pi.f * N.b.NR)
  
  #========================================================================
  
  #7
  mu.a.R <- 0.1 + log(2) * mu.l.a + 0.1 * mu.l.a.R
  mu.a.NR.f <- 0.1 + log(2) * mu.l.a + 0.5 * mu.l.a.NR
  mu.b.R <- 0.1 + log(2) * mu.l.b + 0.1 * mu.l.b.R
  mu.b.NR.f <- 0.1 + log(2) * mu.l.b + 0.5 * mu.l.b.NR
  
  #========================================================================
  
  #8
  mu.d1 <- gam.a * mu.a.R + (1-gam.a)*mu.a.NR.f
  mu.d3 <- gam.b * mu.b.R + (1-gam.b)*mu.b.NR.f
  
  del1 <- mu.d1
  del2 <- mu.d1 - mu.d3
  
  #========================================================================
  
  #9
  sigma.d1 <- runif(1, 1, 10)
  sigma.d3 <- sqrt(r.1 * sigma.d1^2)
  
  nu.d1 <- runif(1, 1, 10)
  nu.d3 <- r.2 * nu.d1
  
  Y.a.R <- as.numeric(rsn(N.a.R, mu.a.R, sigma.d1, nu.d1))
  Y.a.NR.f <- as.numeric(rsn(N.a.NR.f, mu.a.NR.f, sigma.d1, nu.d1))
  Y.b.R <- as.numeric(rsn(N.b.R, mu.b.R, sigma.d1, nu.d1))
  Y.b.NR.f <- as.numeric(rsn(N.b.NR.f, mu.b.NR.f, sigma.d1, nu.d1))
  
  #========================================================================
  
  #10
  
  #need to check W.d1 and W.d3
  
  w.d1 <- (4 - 2*gam.a) 
  Y.bar.d1 <- w.d1 * (sum(Y.a.R) + sum(Y.a.NR.f)) / (w.d1*(N.a.R + N.a.NR.f))
  var.cap.Y.bar.d1 <- (sum((Y.a.R - Y.bar.d1)^2) + sum((Y.a.NR.f - Y.bar.d1)^2)) * w.d1^2 / (N.a.R + N.a.NR.f)^2 
  
  w.d3 <- (4 - 2*gam.b) 
  Y.bar.d3 <- w.d3*(sum(Y.b.R) + sum(Y.b.NR.f))/(w.d3 * (N.b.R + N.b.NR.f))
  var.cap.Y.bar.d3 <- (sum((Y.b.R - Y.bar.d3)^2) + sum((Y.b.NR.f - Y.bar.d3)^2)) * w.d3^2 / (N.b.R + N.b.NR.f)^2
  
  #========================================================================
  
  #11
  del1.cap <- Y.bar.d1
  w1.cap <- qnorm(1 - alpha / 2) * sqrt(var.cap.Y.bar.d1)
  
  del2.cap <- Y.bar.d1 - Y.bar.d3
  w2.cap <- qnorm(1 - alpha /2) * sqrt(var.cap.Y.bar.d1 + var.cap.Y.bar.d3)
  
  CP.del1.a[i] <- (del1 <= del1.cap + w1.cap) & (del1 >= del1.cap - w1.cap)
  CP.del2.a[i] <- (del2 <= del2.cap + w2.cap) & (del2 >= del2.cap - w2.cap)
  psi.cap.del1.a[i] <- w1.cap / sqrt(var.cap.Y.bar.d1)
  psi.cap.del2.a[i] <- w2.cap / sqrt(var.cap.Y.bar.d1 + var.cap.Y.bar.d3)
  
  print(paste0("We are at step ", i))
}

CP.del1 <- mean(CP.del1.a)
CP.del2 <- mean(CP.del2.a)

psi.cap.del1 <- mean(psi.cap.del1.a)
psi.cap.del2 <- mean(psi.cap.del2.a)

