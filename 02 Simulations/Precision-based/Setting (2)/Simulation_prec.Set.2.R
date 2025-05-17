library(sn)                   #for skew normal
library(fGarch)               #for skew-normal parameter estimation

######=================================================================

rm(list = ls())

###------------------------------
###calculating sample sizes

samp.size.prec.s2 <- function(gam.11, gam.12, psi.d1, nu.d1, r1 = 1, r2 = 1, alpha = 0.05)
{
  xi.d1 <- nu.d1 / sqrt(1+nu.d1^2)
  kap.d1 <- (1 - 2*xi.d1^2 / pi) 
  
  ret <- ( (4 - 2*gam.11) + r1 * r2 * (4 - 2*gam.12) ) * qnorm(1-alpha/2)^2 / (psi.d1^2 * kap.d1) 
  
  return(ceiling(ret))
}

###------------------------------
###drawing random samples from truncated skew-normal distribution

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

nu.2 <- function(nu.1, r.3){
  
  sign.nu <- sign(nu.1)
  A <- (1 - r.3) * pi / 2 + r.3 * nu.1^2 / (1+nu.1^2)
  nu.2 <- sign.nu * sqrt(abs(A / (1 - A)))
  
  return(nu.2)
}

#=============================================================================
#the main function

sim.fun.prec.set.2 <- function(psi, gam.11, gam.12, nu.d1, r1 = 1, r2 = 1,
                         rep = 3e3, N = NA, alpha = 0.05, pi.11 = 0.5, pi.22 = 0.5)
{
  
  count <- 0                              
  CP.del.a <- numeric()                    #.a denotes 'all'               
  psi.cap.del.a <- numeric()               #.a denotes 'all'
  
  #1------
  
  if (is.na(N)){N = samp.size.prec.s2(gam.11, gam.12, psi, nu.d1)}
  
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
    if (N.1 == 0) {next}
    
    L.1 <- as.numeric(rsn(N.1, mu.l.1, sigma.l.1, nu.l.1))
    
    # ~ ~ ~ ~ ~
    
    N.2 <- ceiling(N * pi.22)
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
    
    # ~ ~ ~ ~ ~
    
    m.2.R <- r.trun.sn(N.2.R, mu.l.2, sigma.l.2, nu.l.2, ita, Inf)
    m.2.NR <- r.trun.sn(N.2.NR, mu.l.2, sigma.l.2, nu.l.2, -Inf, ita)
    
    if (length(m.2.R)  != N.2.R  | length(m.2.NR)  != N.2.NR) {next}
    
    foo <- snormFit(m.2.R)$par
    mu.m.2.R <- foo[[1]]
    
    foo <- snormFit(m.2.NR)$par
    mu.m.2.NR <- foo[[1]]
    
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
    
    mu.2.R.1 <- 0.1 + log(2) * mu.l.2 + 0.1 * mu.m.2.R
    mu.2.NR.2 <- 0.1 + log(2) * mu.l.2 + 0.5 * mu.m.2.NR
    
    #8------
    mu.d1 <- gam.11 * mu.1.R.1 + (1-gam.11) * mu.1.NR.2
    
    # ~ ~ ~ ~ ~
    
    mu.d3 <- gam.12 * mu.2.R.1 + (1-gam.12) * mu.2.NR.2
    
    del <- mu.d1 - mu.d3
    
    #9------
    sigma.d1 <- runif(1, 3, 6)
    
    Y.1.R.1 <- as.numeric(rsn(N.1.R, mu.1.R.1, sigma.d1, nu.d1))
    Y.1.NR.2 <- as.numeric(rsn(N.1.NR.2, mu.1.NR.2, sigma.d1, nu.d1))
    
    # ~ ~ ~ ~ ~
    
    sigma.d3 <- sqrt(r1) * sigma.d1
    nu.d3 <- nu.2(nu.d1, r2)
    
    Y.2.R.1 <- as.numeric(rsn(N.2.R, mu.2.R.1, sigma.d3, nu.d3))
    Y.2.NR.2 <- as.numeric(rsn(N.2.NR.2, mu.2.NR.2, sigma.d3, nu.d3))
    
    #10------
    Y.bar.d1 <- (sum(2* Y.1.R.1) + sum(4 * Y.1.NR.2)) / (2*length(Y.1.R.1) + 4*length(Y.1.NR.2))
    var.cap.Y.bar.d1 <- (   sum(4 * (Y.1.R.1 - Y.bar.d1)^2) + sum(16 * (Y.1.NR.2 - Y.bar.d1)^2)   ) / N^2
    
    # ~ ~ ~ ~ ~
    
    Y.bar.d3 <- (sum(2* Y.2.R.1) + sum(4 * Y.2.NR.2)) / (2*length(Y.2.R.1) + 4*length(Y.2.NR.2))
    var.cap.Y.bar.d3 <- (   sum(4 * (Y.2.R.1 - Y.bar.d3)^2) + sum(16 * (Y.2.NR.2 - Y.bar.d3)^2)   ) / N^2
    
    #11------
    del.cap <- (Y.bar.d1 - Y.bar.d3) - ( sigma.d1 * ( nu.d1 / sqrt(1 + nu.d1^2) ) - sigma.d3 * ( nu.d3 / sqrt(1 + nu.d3^2) ) ) * sqrt(2/pi)
    w.cap <- qnorm(1 - alpha / 2) * sqrt(var.cap.Y.bar.d1 + var.cap.Y.bar.d3)
    
    sig.cap.1.R.1 <- snormFit(Y.1.R.1)[[1]][2]
    w.1.R.1 <- N.1.R / (N.1.R + N.1.NR.2)
    
    sig.cap.1.NR.2 <- snormFit(Y.1.NR.2)[[1]][2]
    w.1.NR.2 <- N.1.NR.2 / (N.1.R + N.1.NR.2)
    
    sig.d1.cap <- sig.cap.1.R.1 * w.1.R.1 + sig.cap.1.NR.2 * w.1.NR.2
    sd.d1.cap <- sig.d1.cap * sqrt(1 - 2 * nu.d1^2 / (pi * (1+nu.d1^2))  )
    
    #12------
    CP.del <- (del <= del.cap + w.cap) & (del >= del.cap - w.cap)
    psi.cap.del <- w.cap / sd.d1.cap
    
    #---------------------
    #store
    CP.del.a <- c(CP.del.a, CP.del)
    psi.cap.del.a <- c(psi.cap.del.a, psi.cap.del)
    
    #---------------------
    #aesthetics
    if ((i %% 100 == 0) & (i < rep))
    {
      print(paste0("We just finished ", i, " steps."))
    }
    
    count <- count + 1
  }
  
  print(paste0("We just finished ", rep/1000, "k steps."))
  
  ret <- list(Coverage_Prob = mean(na.omit(CP.del.a)), pis_hat = mean(na.omit(psi.cap.del.a)),
              Sample_Size = N, Effective_loops = count)
  
  return(ret)
}

#===============================================================================
#tabulating values

# psi.vec <- c(0.3, 0.4, 0.5)
# gam.11.vec <- c(0.35, 0.5, 0.65)
# gam.12.vec <- c(0.35, 0.5, 0.65)
# nu.d1.vec <- c(-25, -10, -5, -1, -0.5, -0.25, 0, 0.25, 0.5, 1, 5, 10, 25)
# 
# store.cov_prob <- array(dim = c(length(psi.vec),length(gam.11.vec),length(nu.d1.vec)))
# store.psi_hat <- array(dim = c(length(psi.vec),length(gam.11.vec),length(nu.d1.vec)))
# store.n <- array(dim = c(length(psi.vec),length(gam.11.vec),length(nu.d1.vec)))
# 
# for (i.psi in 1:length(psi.vec))
# {
#   for (i.gam in 1:length(gam.11.vec))
#   {
#     for (i.nu.d1 in 1:length(nu.d1.vec))
#     {
#       foo <- sim.fun.prec.set.2(psi.vec[i.psi], gam.11.vec[i.gam],
#                           gam.12.vec[i.gam], nu.d1.vec[i.nu.d1], rep = 1500)
#       
#       store.cov_prob[i.psi, i.gam, i.nu.d1] <- foo$Coverage_Prob
#       store.psi_hat[i.psi, i.gam, i.nu.d1] <- foo$pis_hat
#       store.n[i.psi, i.gam, i.nu.d1] <- foo$Sample_Size
#       
#       #aesthetics
#       print("--------------------------------")
#       print(paste0("We are at psi = ", psi.vec[i.psi], ", gam.a = ", gam.11.vec[i.gam], ", and nu.d1 = ", nu.d1.vec[i.nu.d1],"."))
#       print("--------------------------------")
#     }
#   }
# }
# 
# save(store.cov_prob, store.psi_hat, store.n, file = "Simulation_prec.Set.2.Rdata")

#-------------------------------------------------------------------------------
# Using parallel core

library(foreach)
library(doParallel)

psi.vec <- c(0.3, 0.4, 0.5)
gam.11.vec <- c(0.35, 0.5, 0.65)
gam.12.vec <- c(0.35, 0.5, 0.65)
nu.d1.vec <- c(-25, -10, -5, -1, -0.5, -0.25, 0, 0.25, 0.5, 1, 5, 10, 25)

num_cores <- 13

cl <- makeCluster(num_cores)
registerDoParallel(cl)

parallel_function<- function(i.nu.d1) {
  
  library(sn)
  library(fGarch)
  
  mat.cur.cov_prob <- matrix(0, ncol = length(gam.11.vec), nrow = length(psi.vec))
  mat.cur.psi_hat <- matrix(0, ncol = length(gam.11.vec), nrow = length(psi.vec))
  mat.cur.n <- matrix(0, ncol = length(gam.11.vec), nrow = length(psi.vec))
  
  for (i.psi in 1:length(psi.vec))
  {
    for (i.gam in 1:length(gam.11.vec))
    {
      foo <- sim.fun.prec.set.2(psi.vec[i.psi], gam.11.vec[i.gam],
                          gam.12.vec[i.gam], nu.d1.vec[i.nu.d1], rep = 1500)
      
      mat.cur.cov_prob[i.psi, i.gam] <- foo$Coverage_Prob
      mat.cur.psi_hat[i.psi, i.gam] <- foo$pis_hat
      mat.cur.n[i.psi, i.gam] <- foo$Sample_Size
      
    }
  }
  
  list.out <- list(mat.cur.cov_prob, mat.cur.psi_hat, mat.cur.n)
  return(list.out)
  
}

list.nu.d1 <- foreach(i.nu.d1 = 1:length(nu.d1.vec)) %dopar% {
  parallel_function(i.nu.d1)
}

stopCluster(cl)

store.n <- store.psi_hat <- store.cov_prob <- array(dim = c(length(psi.vec), length(gam.11.vec), length(nu.d1.vec)))

for (i in 1:length(nu.d1.vec)){
  store.cov_prob[,,i] <- list.nu.d1[[i]][[1]]
  store.psi_hat[,,i] <- list.nu.d1[[i]][[2]]
  store.n[,,i] <- list.nu.d1[[i]][[3]]
}

save(store.cov_prob, store.psi_hat, store.n, file = "Simulation_prec.Set.2.Rdata")

####================================================================

rm(list = ls())
load("Simulation_prec.Set.2.Rdata")

#-- -- -- -- -- 

#coverage probability
plot(density(as.vector(store.cov_prob)),
     main = "Density of the coverage pribabilities",
     ylab = "Density",
     xlab = expression(alpha))

#-- -- -- -- --      

#sample_size:
#(Done in a different R file)

#-- -- -- -- -- 
#psi.hat
store.psi_hat
#(other visualizations are in different R file)

####================================================================
# For tabulation in latex

n.nu <- dim(store.cov_prob)[3]

mat.latex <- matrix(0, ncol = 28, nrow = n.nu)
mat.latex[,1] <- c(-25, -10, -5, -1, -0.5, -0.25, 0, 0.25, 0.5, 1, 5, 10, 25)

for (i in 1:n.nu) {
  foo.mat <- matrix(0, ncol = 9, nrow = 3)
  
  foo.mat[1,] <- as.vector(t(store.n[,,i]))
  foo.mat[2,] <- as.vector(t(store.cov_prob[,,i]))
  foo.mat[3,] <- as.vector(t(store.psi_hat[,,i]))
  
  mat.latex[i,2:28] <- as.vector(foo.mat)
}

mat.latex <- round(mat.latex, 3)

df.latex <- as.data.frame(mat.latex)

df.upper <- df.latex[,1:13]
df.lower <- df.latex[,14:28]

# Now use chatgpt to get the latex code
