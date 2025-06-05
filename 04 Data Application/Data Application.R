rm(list = ls())

library(fGarch)

data <- read.csv("obj2_weight_difference.csv")
wgt <- na.omit(data$weight_difference_kg)

nu <- snormFit(wgt)$par[3]

## ============================================================================
## Precision Based Setting (1)

samp.size.prec.s1 <- function(gam.11, psi.d1, nu.d1, alpha = 0.05)
{
  xi.d1 <- nu.d1 / sqrt(1+nu.d1^2)
  kap.d1 <- (1 - 2*xi.d1^2 / pi) 
  ret <- (4 - 2*gam.11) * qnorm(1-alpha/2)^2 / (psi.d1^2 * kap.d1) 
  
  return(ceiling(ret))
}

psi.vec <- c(0.3, 0.4, 0.5)
gam.11.vec <- c(0.35, 0.5, 0.65)

store.n.prec.s1 <- matrix(0, nrow = length(psi.vec), ncol = length(gam.11.vec))

for (i.psi in 1:length(psi.vec))
{
  for (i.gam in 1:length(gam.11.vec))
  {
    gam.11 <- gam.11.vec[i.gam]
    psi <- psi.vec[i.psi]
    
    #-- - -- --
    
    store.n.prec.s1[i.psi, i.gam] <- samp.size.prec.s1(gam.11, psi, nu)
  }
}

store.n.prec.s1 <- as.data.frame(store.n.prec.s1)
colnames(store.n.prec.s1) <- c("gam.35", "gam.5", "gam.65")
rownames(store.n.prec.s1) <- c("psi.3", "psi.4", "psi.5")
store.n.prec.s1

## ============================================================================
## Precision Based Setting (2)

samp.size.prec.s2 <- function(gam.11, gam.12, psi.d1, nu.d1, r1 = 1, r2 = 1, alpha = 0.05)
{
  xi.d1 <- nu.d1 / sqrt(1+nu.d1^2)
  kap.d1 <- (1 - 2*xi.d1^2 / pi) 
  
  ret <- ( (4 - 2*gam.11) + r1 * r2 * (4 - 2*gam.12) ) * qnorm(1-alpha/2)^2 / (psi.d1^2 * kap.d1) 
  
  return(ceiling(ret))
}

psi.vec <- c(0.3, 0.4, 0.5)
gam.vec <- c(0.35, 0.5, 0.65)

store.n.prec.s2 <- matrix(0, nrow = length(psi.vec), ncol = length(gam.vec))

for (i.psi in 1:length(psi.vec))
{
  for (i.gam in 1:length(gam.11.vec))
  {
    gam.11 <- gam.vec[i.gam]
    gam.12 <- gam.vec[i.gam]
    psi <- psi.vec[i.psi]
    
    #-- - -- --
    
    store.n.prec.s2[i.psi, i.gam] <- samp.size.prec.s2(gam.11, gam.12, psi, nu)
  }
}

store.n.prec.s2 <- as.data.frame(store.n.prec.s2)
colnames(store.n.prec.s2) <- c("gam.35", "gam.5", "gam.65")
rownames(store.n.prec.s2) <- c("psi.3", "psi.4", "psi.5")
store.n.prec.s2

## ============================================================================
## Power Based

samp.size.pow.s2 <- function(gam.11 = 0.35, gam.12 = 0.35, nu.d1 = 1, eff.size = 0.2,
                             pow = 0.8, alpha = 0.05, r1 = 1, r2 = 1)
{
  xi.d1 <- nu.d1 / sqrt(1+nu.d1^2)
  kap.d1 <- (1 - 2*xi.d1^2 / pi) 
  ret <- 2 * (qnorm(pow) + qnorm(1- alpha/2))^2 * ((4 - 2*gam.11) + r1 * r2 * (4 - 2*gam.12)) * kap.d1 / (eff.size^2 * (1+r1))
  return(ceiling(ret))
}

eff.vec <- c(0.15, 0.2, 0.25)
gam.vec <- c(0.35, 0.5, 0.65)

store.n.pow <- matrix(0, nrow = length(eff.vec), ncol = length(gam.vec))

for (i.eff in 1:length(eff.vec))
{
  for (i.gam in 1:length(gam.11.vec))
  {
    gam.11 <- gam.vec[i.gam]
    gam.12 <- gam.vec[i.gam]
    eff <- eff.vec[i.eff]
    
    #-- - -- --
    
    store.n.pow[i.eff, i.gam] <- samp.size.pow.s2(gam.11, gam.12, nu, eff)
  }
}

store.n.pow <- as.data.frame(store.n.pow)
colnames(store.n.pow) <- c("gam.35", "gam.5", "gam.65")
rownames(store.n.pow) <- c("eff.15", "eff.2", "eff.25")

## ============================================================================
## Power Based

store.n.prec.s1
store.n.prec.s2
store.n.pow
